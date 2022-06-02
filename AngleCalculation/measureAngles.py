import prody
from datetime import datetime
import argparse
import numpy

parser = argparse.ArgumentParser(description='Program to calculate the phi, psi or both angles in a selection or '
                                             'set of selections.')
parser.add_argument('-i', '--in_file', type=str, required=True, help='Path, either absolute or relative, to the input '
                                                                     'file')
# parser.add_argument('-a', '--angles', type=str, required=False, help='"phi", "psi" or "both" (without ""). Default: '
#                                                                      'both')

arg = parser.parse_args()

dcdName = []
outName = 'outFile.out'

print('Reading variables from input file')
inFile = open(arg.in_file, 'r')

print('Reading input file. Stand by')

for line in inFile:
    l = line.strip().split()
    if 'dcd' in l[0]:
        dcdName.append(l[1])
    elif 'sel' in l[0]:
        selName = ' '.join(l[1:])
    elif 'ref' in l[0]:
        if len(l[1:]) > 1:
            refSel = ' '.join(l[1:])
        else:
            refSel = l[1]
    elif 'pdb' in l[0]:
        pdbName = l[1]
    elif 'out' in l[0]:
        outName = l[1]

inFile.close()

# Ask whether the user forgot a keyword
if 'protein' not in selName and 'chain' not in selName:
    answer = str(input('Did you mean "protein/chain X and ' + selName + '"? If so, please specify. Otherwise, type "NO"'))
    if answer != 'NO' and answer != '"NO"':
        if 'and' not in answer:
            selName = answer + ' and ' + selName
        else:
            selName = answer + ' ' + selName
        if 'resid' not in answer:
            answer = answer + 'and resid'
    else:
        answer = selName.split('resid')[0]
        if len(answer) < 5:
            answer = 'resid'

pdb = prody.parsePDB(pdbName)
traj = prody.Trajectory(dcdName[0])

print('Starting calculations and setting up stuff. You\´re free to go make some coffee or something')
t1 = datetime.now()

# Check for unique residues and split the measurement
unis = numpy.unique(pdb.select(selName).getResnum())
if len(unis) > 1:
    residues = []
    for res in unis:
        residues.append(answer + str(res))
else:
    residues = [selName]

# If there are more dcd files, add them to the trajectory
if len(dcdName) > 1:
    for dcd in dcdName[1:]:
        traj.addFile(dcd)

traj.link(pdb)
traj.setCoords(pdb)
traj.setAtoms(pdb.select(refSel))
dih = numpy.zeros((len(residues) + 1, len(traj))) # Array where angles are to be stored
for i, frame in enumerate(traj):
    frame.superpose()
    for j, res in enumerate(residues):
        if j != 0:
            dih[j+1, i] = prody.calcPsi(res)
        elif j == 0:
            dih[j+1, i] = prody.calcPsi(res)
            dih[j, i] = prody.calcPhi(res)

# Cast as string to be able to write to file easier
dih = dih.astype('str')
unisDummy = unis.copy().astype('str')

outFile = open(outName, 'w+')
outFile.write('#Frame \t {0} \n'.format('\t '.join(unisDummy)))
for k, vals in enumerate(dih):
    outFile.write('{0:7} \t {1} \n'.format(k, '\t '.join(vals)))
outFile.close()

t2 = datetime.now()
print('Done. Completion took {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
