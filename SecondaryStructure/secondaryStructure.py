import mdtraj
from prody import parsePDB
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = 'Calculate the collective dipole of the given selection.')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')

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
dssp = mdtraj.compute_dssp(traj)

with open(outName, mode='w') as f1:
    for line in dssp:
        f1.write(','.join(line.astype(str)) + '\n')

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))