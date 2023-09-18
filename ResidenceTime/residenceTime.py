
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
                spaceBins = int(l[1])
            elif l[0].lower() == 'rad':
                rad = float(l[1])
            elif l[0].lower() == 'thr':
                timeBins = l[1:] # Can be a space-separated list. Each element is a number of minimum frames that define a "time-wise bin"
                            # E.g. If 5, 100, 500 is given, data about residence time in intervals < 5 | 5 - 100 | 100 - 500 | 500 > will be returned.
        else:
            print("--bins (-b) flag was either not set or was set to False. Data won't be binned.")
    else:
        pass

inFile.close()

pdb = prody.parsePDB(pdbName)
dcd = prody.Trajectory(dcdName[0])
if len(dcdName) > 1:
    for d in dcdName[1:]:
        dcd.addFile(d)

dcd.link(pdb)
dcd.setCoords(pdb)
dcd.setAtoms(pdb.select(refName)) # refName = Selection used when aligning frames (frame.superpose())
atomsSelection = pdb.select(selName)

moleculeCharacterization = numpy.zeros((len(atomsSelection), len(timeBins), spaceBins))