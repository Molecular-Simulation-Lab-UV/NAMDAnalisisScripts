import prody
import socket
from datetime import datetime
import argparse
import numpy

angle = 'both'

parser = argparse.ArgumentParser(description='Program to calculate the phi, psi or both angles in a selection or '
                                             'set of selections.')
parser.add_argument('-i', '--in_file', type=str, required=True, help='Path, either absolute or relative, to the input '
                                                                     'file')
parser.add_argument('-a', '--angles', type=str, required=False, help='"phi", "psi" or "both" (without ""). Default: '
                                                                      'both')

arg = parser.parse_args()

dcdName = []
outName = 'outFile.out'

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
            if len(l[1:]) > 1:
                refSel = ' '.join(l[1:])
            else:
                refSel = l[1]
        elif 'chain' in l[0]:
            chainID = l[1]
        elif 'pdb' in l[0]:
            pdbName = l[1]
        elif 'out' in l[0]:
            outName = l[1]
    except:
        pass

inFile.close()

pdb = prody.parsePDB(pdbName)
traj = prody.Trajectory(dcdName[0])

if 'resid ' not in selName:
    selName = 'protein and chain ' + chainID + ' and resid ' + selName
elif 'resid' in selName:
    selName = 'protein and chain ' + chainID + ' and ' + selName

print('Starting calculations and setting up stuff. You\'re free to go make some coffee or something')
t1 = datetime.now()

# Check for unique residues and split the measurement
selUniverse = pdb.getHierView()
unis = numpy.unique(pdb.select(selName).getResnums())

# If there are more dcd files, add them to the trajectory
if len(dcdName) > 1:
    for dcd in dcdName[1:]:
        traj.addFile(dcd)

traj.link(pdb)
traj.setCoords(pdb)
traj.setAtoms(pdb.select(refSel))

#Arrays for storing the calculated values

if len(traj) > 1e4 and len(unis) > 3 and angle == 'both' and socket.gethostname != 'don-elias':
    print('Better run this on the server. Otherwise memory will probably go (╯°□°)╯︵ ┻━┻ \n')
    print('Recommendation is to Ctrl + C this process and run it on don-elias.')

dih = numpy.zeros((len(traj), len(unis))) 
if angle == 'both':
    dih2 = numpy.zeros((len(traj), len(unis)))

if angle == 'phi':
    for i, frame in enumerate(traj):
        frame.superpose()
        for j, res in enumerate(unis):
            dih[i, j] = prody.calcPhi(selUniverse[chainID,chainID,res])
elif angle == 'psi':
    for i, frame in enumerate(traj):
        frame.superpose()
        for j, res in enumerate(unis):
            dih[i, j] = prody.calcPsi(selUniverse[chainID,chainID,res])
else:
    for i, frame in enumerate(traj):
        frame.superpose()
        for j, res in enumerate(unis):
            dih[i, j] = prody.calcPhi(selUniverse[chainID,chainID,res])
            dih2[i, j] = prody.calcPsi(selUniverse[chainID,chainID,res])

# Cast as string to be able to write to file easier
dih = dih.astype('str')
if angle == 'both':
    dih2 = dih2.astype('str')

unis = unis.copy().astype('str')

if '.out' in outName:
    outName = outName.split('.out')[0]
if angle != 'both':
    outFile = open(outName + '.' + angle + '.out', 'w+')
    outFile.write('#Frame \t resids {0} \n'.format('\t '.join(unis)))
    for k, vals in enumerate(dih):
        outFile.write('{0:7} \t {1} \n'.format(k, '\t '.join(vals)))
    outFile.close()
else:
    # Write phi bond angles
    outFile1 = open(outName + '.phi.out', 'w+')
    outFile1.write('#Frame \t resids {0} \n'.format('\t '.join(unis)))
    for k, vals in enumerate(dih):
        outFile1.write('{0:7} \t {1} \n'.format(k, '\t '.join(vals)))
    outFile1.close()
    # Write psi bond angles
    outFile2 = open(outName + '.psi.out', 'w+')
    outFile2.write('#Frame \t resids {0} \n'.format('\t '.join(unis)))
    for k, vals in enumerate(dih2):
        outFile2.write('{0:7} \t {1} \n'.format(k, '\t '.join(vals)))
    outFile2.close()


t2 = datetime.now()
print('Done (¬‿¬)=b. Completion took {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
