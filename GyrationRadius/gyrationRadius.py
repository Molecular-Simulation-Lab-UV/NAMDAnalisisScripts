import prody
import numpy as np
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = "Calculate the permeation events of selection's atoms through a pore.")
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

dcdName = []
outName = 'outFile.out'

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
    if 'dcd' in l[0].lower(): # Any amount of .dcd files. Multiple lines, in order up-down, can be input. Only requisite is "dcd" should appear as keyword.
        dcdName.append(l[1])
    elif l[0].lower() == 'pdb':
        pdbName = l[1]
    elif l[0].lower() == 'ref': # Alignment selection for superposition of the .dcd frames onto the pdb conformation.
        if len(l[1:]) > 1:
            refName = ' '.join(l[1:])
        else:
            refName = l[1]
    elif 'sel' in l[0].lower(): # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if len(l[1:]) > 1:
            selName = ' '.join(l[1:])
        else:
            selName = l[1]
    elif l[0].lower() == 'out': # Path to the output file, name included.
        outName = l[1]

inFile.close()

pdb = prody.parsePDB(pdbName)
traj = prody.Trajectory(dcdName[0])
if len(dcdName) > 1:
    for d in dcdName[1:]:
        traj.addFile(d)

t1 = datetime.now()

traj.link(pdb)
traj.setCoords(pdb)
traj.setAtoms(pdb.select(refName)) # refName = Selection used when aligning frames (frame.superpose())
sel = pdb.select(selName)

gyration = np.zeros(len(traj))

for f, frame in enumerate(traj):
    frame.superpose()
    gyration[f] = prody.calcGyradius(sel)

outFile = open(outName, 'w+')
outFile.write('# Frame \t |\t Gyration Radius \n')

for f, vals in enumerate(gyration):
    outFile.write('{0} \t\t {1} \n'.format(f, vals))
    
outFile.close()

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))