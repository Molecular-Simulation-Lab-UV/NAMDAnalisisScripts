"""
This script requires a psf for the partial charges. For this, a matching pair of files, i.e. a .pdb and a .psf, are required.
This script was made in 2022; please consider that in case of eventual updates/changes to file formats.
We take the partial charges for the designated atoms in the input file, and proceed to calculate the dipole moment of that selection.
WARNING: ONLY WORKS FOR A FIXED SELECTION! SOMETHING LIKE "water and z > 4 and z < 10" WON'T WORK IF UPDATED AT LEAST ONCE.
"""

import prody
import numpy
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = 'Calculate the collective dipole of the given selection.')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')
parser.add_argument('-a', '--average', action = 'store_true', required = False, default = False, help = 'Calculate the average for the trajectory if included (with "-a" option), frame-by-frame-wise if not included. Default: False')
parser.add_argument('-b', '--bins', action = 'store_true', required = False, default = False, help = 'Return dipole of selection "sel" binned across a cylinder. Default: False')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

dcdName = []
outName = 'outFile.out'

avg = arg.average


# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
    if len(l) > 0:
        if 'dcd' in l[0]: # Any amount of .dcd files
            dcdName.append(l[1])
        elif l[0] == 'pdb':
            pdbName = l[1]
        elif l[0] == 'psf':
            psfName = l[1]
        elif l[0] == 'ref': # Alignment selection
            if len(l[1:]) > 1:
                refName = ' '.join(l[1:])
            else:
                refName = l[1]
        elif 'sel' in l[0]: # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2"
            if len(l[1:]) > 1:
                selName = ' '.join(l[1:])
            else:
                selName = l[1]
        elif l[0] == 'out':
            outName = l[1]
        if arg.bins == True:
            if l[0].lower() == 'zmin':
                zMin = float(l[1])
            elif l[0].lower() == 'zmax':
                zMax = float(l[1])
            elif l[0].lower() == 'nbins':
                nBins = int(l[1])
            elif l[0].lower() == 'rad':
                rad = float(l[1])
        else:
            print("--bins (-b) flag was either not set or was set to False. Data won't be binned.")
    else:
        pass

inFile.close()

# Setting up structure and linking the trajectory files to it.
pdb = prody.parsePDB(pdbName) # Import and parse the .pdb
atomSystem = prody.parsePSF(psfName, ag=pdb) # Import the .psf and associate the .pdb data to it
dcd = prody.Trajectory(dcdName[0])
if len(dcdName) > 1: # Check for and add corresponding dcd files to the trajectory (concatenate them)
    for dcd1 in dcdName[1:]:
        dcd.addFile(dcd1)

atoms = atomSystem.select(selName).getIndices()
if len(atoms) < 2:
    print('************************')
    print('Something went wrong. You need at least 2 atoms to calculate dipoles!!! (╯°□°)╯︵ ┻━┻')
    print('************************')
    exit()

dcd.link(pdb)
dcd.setCoords(pdb)
dcd.setAtoms(pdb.select(refName))

t1 = datetime.now()
# Calculate and return dipole

if not arg.bins:
    
    sel = atomSystem.select(selName)
    charges = sel.getCharges()
    dipoleArray = numpy.zeros(len(dcd))

    for f, frame in enumerate(dcd):
        prody.wrapAtoms(pdb, unitcell = frame.getUnitcell()[:3], center = prody.calcCenter(pdb.select(refName)))
        frame.superpose()
        center = prody.calcCenter(sel, weights=None) # Get the geometrical center of the selection
        coords = sel.getCoords() # Get each of the selection's atoms' position
        dipole = numpy.sum((coords - center)*charges[:,numpy.newaxis], axis=0) # Calculate the dipole as sum(q*(r_i - r_com)); com should be center of mass, but is taken as geometrical center.
        dipoleM = numpy.sqrt(numpy.sum(dipole*dipole)) # Calculate the magnitude as the dot product of the selection's vector with itself, and then apply the square root.
        dipoleArray[f] = dipoleM

    if arg.average:
        dipoleAvg = numpy.average(dipoleArray)
        dipoleStd = numpy.std(dipoleArray)
        print('\nDipole average is {0:6.3f}'.format(dipoleAvg))
        print('Dipole standard deviation is {0:6.3f}'.format(dipoleStd))

    elif not arg.average:
        outFile = open(outName, 'w+')
        outFile.write('#Frame \t Dipole magnitude \n')

        for ff, vals in enumerate(dipoleArray):
            outFile.write('{0} \t {1} \n'.format(ff, vals))
        
        outFile.close()

elif arg.bins:

    L = zMax - zMin
    binSize = L/nBins
    dipoleArray = numpy.zeros((len(dcd), nBins, 3)) # The 3 are the cartesian coordinates of the dipole moment vector
    binArray = numpy.arange(zMin, zMax + binSize, binSize)

    for f, frame in enumerate(dcd):
        prody.wrapAtoms(pdb, unitcell = frame.getUnitcell()[:3], center = prody.calcCenter(pdb.select(refName)))
        frame.superpose()
        sel = pdb.select(f'same residue as ({selName} and (x^2 + y^2) < {rad**2} and z < {zMax} and z > {zMin})')
        charges = sel.getCharges()
        coords = sel.getCoords()
        sliceIndices = numpy.cumsum(numpy.unique(sel.getResindices(), return_counts = True)[1])[:-1]
        slices = numpy.split(coords, sliceIndices)
        geoCenters = numpy.average(slices, axis = 1)
        dipoles = numpy.sum((slices - geoCenters[:,numpy.newaxis,:])*numpy.reshape(charges, (int(len(charges)/sliceIndices[0]), sliceIndices[0], 1)), axis = 1)
        dipoles = dipoles/numpy.linalg.norm(dipoles, axis = 1)[:,numpy.newaxis]
        flags = numpy.argwhere((numpy.sqrt(geoCenters[:,0]**2 + geoCenters[:,1]**2) < rad)) # Cylinder
        zPos = geoCenters[flags,-1]
        inBin = numpy.argwhere((zPos >= binArray[numpy.newaxis, :-1]) & (zPos < binArray[numpy.newaxis, 1:]))

        binnedDipoles = numpy.zeros((len(binArray)-1, 3))
        binUniques, binCounts = numpy.unique(inBin[:,1], return_counts = True)
        diffArray = numpy.flatnonzero(numpy.isin(numpy.arange(0, len(binArray)-1, 1), binUniques, invert = True))
        binCounts = numpy.insert(binCounts, diffArray - numpy.arange(len(diffArray)), 1) # Number 1 doesn't matter, it could be any non-zero value
        numpy.add.at(binnedDipoles, inBin[:,1], dipoles[inBin[:,0]])
        dipoleArray[f] = binnedDipoles/binCounts[:, numpy.newaxis]
        
    if arg.average:
        dipoleAvg = numpy.average(dipoleArray, axis = 0)
        dipoleAvg = dipoleAvg.astype('str')
        dipoleStd = numpy.std(dipoleArray, axis = 0)
        dipoleStd = dipoleStd.astype('str')
        binCenters = (binArray[:-1] + binArray[1:])/2
        binCenters = binCenters.astype('str')
        outFile = open(outName, 'w+')
        outFile.write('# Bin Center (z) \t Dipole average \t Dipole Std\n')
        for ff, vals in enumerate(dipoleAvg):
            outFile.write('{0} \t {1} \t {2} \n'.format(binCenters[ff], '\t'.join(vals), '\t'.join(dipoleStd[ff])))
        
        outFile.close()
    else:
        dipoleArray = dipoleArray.astype('str')
        binCenters = (binArray[:-1] + binArray[1:])/2
        binCenters = binCenters.astype('str')
        outFile = open(outName, 'w+')
        outFile.write('# Frame \t | \t Dipole average per bin, with bin centers (z): \t {}\n'.format(binCenters))
        for ff, vals in enumerate(dipoleArray):
            outFile.write('{0} \t {1} \n'.format(ff, '\t'.join(numpy.concatenate(vals))))
        
        outFile.close()

t2 = datetime.now()
print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
