import mdtraj
from prody import parsePDB
import numpy as np
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='Calculate the secondary structure of the given selection.')
parser.add_argument('-i', '--in_file', type=str, required=True, help='Path, either absolute or relative, to the input file')
parser.add_argument('-a', '--atom', action='store_true', required=False, default=False, help='Calculate atom-wise or residue-wise SASA. Passing the flag changes to per-atom.')

arg = parser.parse_args()

# dcdName = []
outName = 'outFile.out'

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

with open(arg.in_file, 'r') as f:
    for line in f:
        l = line.strip().split()
        if len(l) > 0:
            if 'dcd' in l[0].lower(): # Any amount of .dcd files
                # dcdName.append(l[1])
                dcdName = l[1]
            elif l[0].lower() == 'pdb':
                pdbName = l[1]
            elif l[0].lower() == 'sel': # Alignment selection
                if len(l[1:]) > 1:
                    selName = ' '.join(l[1:])
                else:
                    selName = l[1]
            elif l[0].lower() == 'out':
                outName = l[1]

pdb = parsePDB(pdbName)
sel = pdb.select(selName)
idx = sel.getIndices()

t1 = datetime.now()

traj = mdtraj.load_dcd(dcdName, top=pdbName, atom_indices=idx)
sasa = mdtraj.shrake_rupley(traj=traj, mode='atom' if arg.atom else 'residue')
sasa = np.hstack((sasa, np.sum(sasa, axis=1)[:,np.newaxis]))

with open(outName, mode='w') as f1:
    f1.write('# ' + selName + ' | Last column is total\n')
    for line in sasa:
        f1.write(','.join(line.astype(str)) + '\n')

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))