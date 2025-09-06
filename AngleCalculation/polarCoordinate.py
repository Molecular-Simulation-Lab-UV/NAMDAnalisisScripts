import prody
from Miscellaneous.reader import Reader
from datetime import datetime
import argparse
import numpy as np
from math import pi

parser = argparse.ArgumentParser(description='Program to calculate the phi, psi or both angles in a selection or '
                                             'set of selections.')
parser.add_argument('-i', '--in_file', type=str, required=True, help='Path, either absolute or relative, to the input '
                                                                     'file')

arg = parser.parse_args()
r = Reader(arg.in_file)

pdb = prody.parsePDB(r.variables['pdb'])
traj = prody.Trajectory('dcd')
for dcd in r.variables['dcd']:
        traj.addFile(dcd)
traj.link(pdb)
traj.setCoords(pdb)
traj.setAtoms(pdb.select(r.variables['ref']))

print('Starting calculations and setting stuff up. You\'re free to go make some coffee or something')

t1 = datetime.now()

# TODO: Sel should be conditional. COM or individual atoms (e.g. resid 1), a choice must be made
sel = pdb.select(r.variables['sel'])
theta = np.zeros((len(traj), len(sel)))

for f, frame in enumerate(traj):
    frame.superpose()
    pos = sel.getCoords()
    norm_pos = np.linalg.norm(pos, axis=1)
    projection_x = np.divide(np.dot(pos, [1, 0, 0]), norm_pos)
    projection_y = np.divide(np.dot(pos, [0, 1, 0]), norm_pos)
    theta[f] = 180*np.arctan2(pos[:,1], pos[:,0])/pi

with open(r.variables['out'], 'w+') as out_file:
    out_file.write('# Frame \t | \t {0} \n'.format(r.variables['sel']))
    theta = theta.astype('str')
    for f, vals in enumerate(theta):
        out_file.write('{0:7} \t {1} \n'.format(f, '\t'.join(vals)))
    
t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))