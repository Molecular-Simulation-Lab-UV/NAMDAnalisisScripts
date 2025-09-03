import prody
from datetime import datetime
import argparse
import numpy as np
from math import pi

parser = argparse.ArgumentParser(description='Program to calculate the phi, psi or both angles in a selection or '
                                             'set of selections.')
parser.add_argument('-i', '--in_file', type=str, required=True, help='Path, either absolute or relative, to the input '
                                                                     'file')

arg = parser.parse_args()

dcdName = []
outName = 'outFile.dat'

inFile = open(arg.in_file, 'r')

print('Reading input file. Stand by')

for line in inFile:
    l = line.strip().split()
    try:
        if 'dcd' in l[0]:
            dcdName.append(l[1])
        elif 'sel' in l[0]:
            selName = ' '.join(l[1:])
        elif 'ref' in l[0]:
            refSel = ' '.join(l[1:])
        elif 'pdb' in l[0]:
            pdbName = l[1]
        elif 'out' in l[0]:
            outName = l[1]
    except:
        pass

inFile.close()

pdb = prody.parsePDB(pdbName)
traj = prody.Trajectory(dcdName[0])
if len(dcdName) > 1:
    for dcd in dcdName[1:]:
        traj.addFile(dcd)
traj.link(pdb)
traj.setCoords(pdb)
traj.setAtoms(pdb.select(refSel))

print('Starting calculations and setting stuff up. You\'re free to go make some coffee or something')

t1 = datetime.now()

# TODO: Sel should be conditional. COM or individual atoms (e.g. resid 1), a choice must be made
sel = pdb.select(selName)
theta = np.zeros((len(traj), len(sel)))

for f, frame in enumerate(traj):
    frame.superpose()
    pos = sel.getCoords()
    theta[f] = 180*np.arccos(np.divide(1, np.linalg.norm(pos, axis=1)))/pi

outFile = open(outName, 'w+')
outFile.write('# Frame \t |\t {0} \n'.format(selName))

for f, vals in enumerate(theta):
    outFile.write('{0} \t\t {1} \n'.format(f, vals))
    
outFile.close()

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))