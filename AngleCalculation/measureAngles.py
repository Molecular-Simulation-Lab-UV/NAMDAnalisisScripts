import prody
from datetime import datetime
import argparse
import numpy as np

angle = 'both'

parser = argparse.ArgumentParser(description='Program to calculate the phi, psi or both angles in a selection or '
                                             'set of selections.')
parser.add_argument('-i', '--in_file', type=str, required=True, help='Path, either absolute or relative, to the input '
                                                                     'file')
parser.add_argument('-a', '--angles', type=str, required=False, help='"phi", "psi" or "both" (without ""). Default: '
                                                                      'both')

arg = parser.parse_args()

dcdName = []
outName = 'outFile.dat'

if arg.angles != None:
    angle = arg.angles
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

# Check for unique residues and split the measurement
sel = pdb.select(selName)
view = sel.getHierView()
unis = np.unique(sel.getResnums())
if (np.unique(sel.getSegnames()) == '').all():
    for chain in np.unique(sel.getChids()):
        pdb.select('chain ' + chain).setSegnames('PRO' + chain)
elif (np.unique(sel.getChids()) == '').all():
    for seg in np.unique(sel.getSegnames()):
        pdb.select('segment ' + seg).setChids(seg[-1])

segment = np.unique(sel.getSegnames())[0]
chain = np.unique(sel.getChids())[0]

if view.numSegments() > 1 or view.numChains() > 1:
    print('Please run a single job per segment/chain. Modify your selection and try again')
    exit()

t1 = datetime.now()
sel2 = pdb[segment, chain]
# Arrays for storing the calculated values
dih = np.zeros((len(traj), len(unis))) 
if angle == 'both':
    dih2 = np.zeros((len(traj), len(unis)))

if angle == 'phi':
    for i, frame in enumerate(traj):
        frame.superpose()
        for j, res in enumerate(unis):
            dih[i, j] = prody.calcPhi(sel2[res,], dist=None)
elif angle == 'psi':
    for i, frame in enumerate(traj):
        frame.superpose()
        for j, res in enumerate(unis):
            dih[i, j] = prody.calcPsi(sel2[res,], dist=None)
else:
    for i, frame in enumerate(traj):
        frame.superpose()
        for j, res in enumerate(unis):
            dih[i, j] = prody.calcPhi(sel2[res,], dist=None)
            dih2[i, j] = prody.calcPsi(sel2[res,], dist=None)

# Cast as string to be able to write to file easier
dih = dih.astype('str')
if angle == 'both':
    dih2 = dih2.astype('str')

unis = unis.copy().astype('str')

if '.out' in outName:
	ending = '.out'
elif '.dat' in outName:
	ending = '.dat'
else:
	ending = '.' + outName.split('.')[-1]

outName = outName.split(ending)[0]

if angle != 'both':
    print('\nWriting to output file:\t' + outName + '.' + angle + ending + '\n')
    outFile = open(outName + '.' + angle + ending, 'w+')
    outFile.write('#Frame \t resids {0} \n'.format('\t '.join(unis)))
    for k, vals in enumerate(dih):
        outFile.write('{0:7} \t {1} \n'.format(k, '\t '.join(vals)))
    outFile.close()
else:
    # Write phi bond angles
    print('\nWriting to output file:\t' + outName + '.' + 'phi' + ending)
    outFile1 = open(outName + '.phi' + ending, 'w+')
    outFile1.write('#Frame \t resids {0} \n'.format('\t '.join(unis)))
    for k, vals in enumerate(dih):
        outFile1.write('{0:7} \t {1} \n'.format(k, '\t '.join(vals)))
    outFile1.close()
    # Write psi bond angles
    print('Writing to output file:\t' + outName + '.' + 'psi' + ending + '\n')
    outFile2 = open(outName + '.psi' + ending, 'w+')
    outFile2.write('#Frame \t resids {0} \n'.format('\t '.join(unis)))
    for k, vals in enumerate(dih2):
        outFile2.write('{0:7} \t {1} \n'.format(k, '\t '.join(vals)))
    outFile2.close()


t2 = datetime.now()
print('Done (¬‿¬)=b. Completion took {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
