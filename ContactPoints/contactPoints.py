import prody
import numpy
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = "Calculate the permeation events of selection's atoms through an equivalent cylinder.")
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

thrshld = 5
dcdName = []
outName = 'outFile.out'

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
    if 'dcd' in l[0]: # Any amount of .dcd files
        dcdName.append(l[1])
    elif l[0].lower() == 'pdb':
        pdbName = l[1]
    elif l[0].lower() == 'thr':
        thrshld = l[1]
    elif l[0].lower() == 'ref': # Alignment selection
        if len(l[1:]) > 1:
            refName = ' '.join(l[1:])
        else:
            refName = l[1]
    elif l[0].lower() == 'sel1': # First selection. Contact points for all of the atoms in it will be counted.
        if len(l[1:]) > 1:
            mainSelName = ' '.join(l[1:])
        else:
            mainSelName = l[1]
    elif l[0].lower() == 'sel2': # Second selection. Contact points between sel1 and this one will be counted.
        if len(l[1:]) > 1:
            selName = ' '.join(l[1:])
        else:
            selName = l[1]
    elif l[0].lower() == 'out':
        outName = l[1]

pdb = prody.parsePDB(pdbName)
traj = prody.Trajectory(dcdName[0])
if len(dcdName) > 1:
    for d in dcdName[1:]:
        traj.addFile(d)

traj.link(pdb)
traj.setCoords(pdb)
traj.setAtoms(pdb.select(refName)) # refName = Selection used when aligning frames (frame.superpose())
mainSel = pdb.select(mainSelName)
secondSel = pdb.select(selName)

print('Calculating contact points between selections\n')
print(mainSelName + ' consisting of {0} atoms'.format(len(mainSel)))
print('and')
print(selName + ' consisting of {0} atoms.\n'.format(len(secondSel)))

t1 = datetime.now()
countArray = numpy.zeros((len(traj), 2))
countArray[:,0] = numpy.linspace(1, len(traj), len(traj))

for i, frame in enumerate(traj):
    prody.wrapAtoms(pdb, unitcell = frame.getUnitcell()[:3], center = prody.calcCenter(pdb.select(refName)))
    frame.superpose()
    selContacts = prody.Contacts(mainSel)
    contacts = selContacts.select(thrshld, secondSel)
    # mainSelCoords = mainSel.getCoords()
    # secondSelCoords = secondSel.getCoords()
    # diffArray = mainSelCoords[:, numpy.newaxis, :] - secondSelCoords[numpy.newaxis, :, :]
    # diffArrayNorm = numpy.linalg.norm(diffArray)
    # numContacts = len(numpy.flatnonzero(numpy.argwhere(diffArrayNorm < thrshld)))
    # countArray[i,1] = numContacts
    countArray[i,1] = len(contacts) if len(contacts) > 0 else 0

countArray = countArray.astype('str') # Handy for writing to a file
outFile = open(outName, 'w+')
outFile.write('# Frame \t\t Contacts count\n')
for vals in countArray:
    outFile.write('{0} \n'.format(' \t '.join(vals)))
outFile.close()

t2 = datetime.now()
print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))