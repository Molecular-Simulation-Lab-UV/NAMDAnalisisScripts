"""
This script requires a psf for the partial charges. For this, a matching pair of files, i.e. a .pdb and a .psf, are required.
This script was made in 2022; please consider that in case of eventual updates/changes to file formats.
We take the partial charges for the designated atoms in the input file, and proceed to calculate the dipole moment of that selection.
WARNING: ONLY WORKS FOR A FIXED SELECTION! SOMETHING LIKE "water and z > 4 and z < 10" WON'T WORK IF UPDATED AT LEAST ONCE.
"""

import prody
import numpy
import pandas
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = 'Calculate the collective dipole of the given selection.')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')
parser.add_argument('-a', '--average', type = bool, required = False, help = 'Calculate the average for the trajectory (True), frame-by-frame-wise (False). Default is True (average).')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

dcdName = []
outName = 'outFile.out'

if arg.average != None:
    avg = arg.average
else:
    avg = True

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
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

inFile.close()
# ADD CHECK TO SEE IF THE SELECTION HAS AT LEAST 2 ATOMS!

# Setting up structure and linking the trajectory files to it.

pdb = prody.parsePDB(pdbName)
atomSystem = prody.parsePSF('psfName', ag=pdb)
sel = atomSystem.select(selName)
dcd = prody.Trajectory(dcdName[0])
if len(dcdName) > 1:
    for dcd1 in dcdName[1:]:
        dcd.addFile(dcd1)


atoms = atomSystem.select(selName).getIndices()
if len(atoms) < 2:
    print('************************')
    print('Something went wrong. You need at least 2 atoms to calculate dipoles!!! (╯°□°)╯︵ ┻━┻')
    print('************************')
    exit()

t1 = datetime.now()

charges = sel.getCharges()

dcd.link(pdb)
dcd.setCoords(pdb)
dcd.setAtoms(pdb.select(refName))

# Get atom indexes from the input selection
# Get charges from the dataframe based on the selection
# Calculate and return dipole

dipoleArray = numpy.zeros(len(dcd))

for f, frame in enumerate(dcd):
    frame.superpose()
    center = prody.calcCenter(sel)
    atomPos = sel.getCoords()
    dipole = numpy.sum((atomPos - center)*charges)
    dipole = numpy.sqrt(dipole*dipole)
    dipoleArray[f] = dipole

if arg.average:
    dipoleAvg = numpy.average(dipoleArray)
    dipoleStd = numpy.std(dipoleArray)
    print('\nDipole average is {0:6.3f}'.format(dipoleAvg))
    print('Dipole standard deviation is {0:6.3f}'.format())

# FALTA CALCULAR LA MAGNITUD DE LOS VECTORES EN DIPOLEARRAY

elif not arg.average:
    outFile = open(outName, 'w+')
    outFile.write('#Frame \t Dipole magnitude \n')

    for ff, vals in enumerate(dipoleArray):
        outFile.write('{0} \t {1} \n'.format(ff, vals))
    
    outFile.close()

t2 = datetime.now()
print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
