"""
Script that returns how much water is in a cylinder, split in bins.
"""

import prody
import numpy
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = 'Calculate the collective dipole of the given selection.')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')
parser.add_argument('-a', '--average', action = 'store_true', required = False, default = False, help = 'Calculate the average for the trajectory if included (with "-a" option), frame-by-frame-wise if not included. Default: False')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')
firstFrame = 0

dcdName = []
outName = 'outFile.out'

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
        # elif l[0] == 'psf':
            # psfName = l[1]
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
        if l[0].lower() == 'zmin':
            zMin = float(l[1])
        elif l[0].lower() == 'zmax':
            zMax = float(l[1])
        elif l[0].lower() == 'nbins':
            nBins = int(l[1])
        elif l[0].lower() == 'rad':
            rad = float(l[1])
        elif l[0].lower() == 'ff':
            firstFrame = int(l[1])
    else:
        pass

inFile.close()

# Setting up structure and linking the trajectory files to it.
pdb = prody.parsePDB(pdbName) # Import and parse the .pdb
dcd = prody.Trajectory(dcdName[0])
if len(dcdName) > 1: # Check for and add corresponding dcd files to the trajectory (concatenate them)
    for dcd1 in dcdName[1:]:
        dcd.addFile(dcd1)

dcd.link(pdb)
dcd.setCoords(pdb)
dcd.setAtoms(pdb.select(refName))
dcd.skip(firstFrame)

loadsArray = numpy.zeros((len(dcd), nBins))
binSize = (zMax - zMin)/nBins
binArray = numpy.arange(zMin, zMax + binSize, binSize)

t1 = datetime.now()


for f, frame in enumerate(dcd):
    prody.wrapAtoms(pdb, unitcell = frame.getUnitcell()[:3], center = prody.calcCenter(pdb.select(refName)))
    frame.superpose()
    coords = pdb.select(selName).getCoords()
    flags = numpy.argwhere((numpy.sqrt(coords[:,0]**2 + coords[:,1]**2) < rad) & (coords[:,2] < zMax) & (coords[:,2] > zMin))
    inBin = numpy.argwhere((coords[flags,2] >= binArray[numpy.newaxis, :-1]) & (coords[flags,2] < binArray[numpy.newaxis, 1:]))
    bins, binCounts = numpy.unique(inBin[:,1], return_counts = True)
    loadsArray[f,bins] = binCounts

if arg.average:
    loadsAvg = numpy.average(loadsArray, axis = 0)
    loadsStd = numpy.std(loadsArray, axis = 0)
    binCenters = (binArray[:-1] + binArray[1:])/2
    outFile = open(outName, 'w+')
    outFile.write('# Z\t Loads avg\t Loads std\n')
    for ff, vals in enumerate(loadsAvg):
        outFile.write('{0:2.2f} \t {1:2.5f} \t {2:2.5f} \n'.format(binCenters[ff], vals, loadsStd[ff]))
    
    outFile.close()
else:
    loadsArray = loadsArray.astype('str')
    binCenters = ((binArray[:-1] + binArray[1:])/2).astype('str')
    outFile = open(outName, 'w+')
    outFile.write('# Frame \t {0}\n'.format(binCenters))
    for ff, vals in enumerate(loadsArray):
        outFile.write('{0} \t {1} \n'.format(ff, '\t'.join(vals)))
    
    outFile.close()

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))