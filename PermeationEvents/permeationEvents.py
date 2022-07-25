"""
This program relies on a user-defined cylinder that roughly defines the pore in order to count molecules
going from the inner solvent to the outer solvent throught the pore. To achieve this, the variables
'upperZ', 'lowerZ' and 'rad' are needed in the input file, which differs from other scripts in this repository.
"""

import prody
import numpy
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = 'Calculate the RMSD Matrix of a trajectory for input in the GROMOS++ clustering program.')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')
parser.add_argument('-d', '--direction', type = str, required = False, help = 'Count events in the "up" direction (z+), "down" direction (z-), or "total" events. Default = "total"')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

direction = 'both'

if arg.direction == None or arg.direction.lower() == 'both':
    direction = 'both'
elif arg.direction.lower() == 'up':
    direction = 'up'
elif arg.direction.lower() == 'down':
    direction = 'down'

dcdName = []
outName = 'outFile.out'

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
    if 'dcd' in l[0]: # Any amount of .dcd files
        dcdName.append(l[1])
    elif l[0] == 'pdb':
        pdbName = l[1]
    elif l[0] == 'ref': # Alignment selection
        if len(l[1:]) > 1:
            refName = ' '.join(l[1:])
        else:
            refName = l[1]
    elif l[0].lower() == 'rad':
        rad = float(l[1])
        rad2 = rad**2
    elif l[0].lower() == 'upperz':
        upZ = float(l[1])
    elif l[0].lower() == 'lowerz':
        loZ = float(l[1])
    elif 'sel' in l[0]: # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2"
        if len(l[1:]) > 1:
            selName = ' '.join(l[1:])
        else:
            selName = l[1]
    elif l[0] == 'out':
        outName = l[1]

pdb = prody.parsePDB(pdbName)
traj = prody.DCDFile(dcdName[0])

if len(traj) > 1:
    for d in dcdName[1:]:
        traj.addfile(d)

if direction == 'both':
    NN = 3
elif direction == 'up' or direction == 'down':
    NN = 1

traj.link(pdb)
traj.setCoords(pdb)
traj.setAtoms(pdb.select(refName))
atomsSelection = pdb.select(selName)

# Initialise the dictionary with zeros. Values will change when permeation events are detected, and all three are reset to 0.
# 4 is for extracellular, 2 for inside the cylinder, 1 for intracellular.
dataDict = {}
atomIDs = atomsSelection.getIndices()
for ind in atomIDs:
    dataDict[ind] = 0

flagsOld = numpy.zeros(len(atomsSelection), 2)
permArray = numpy.zeros((len(traj), NN)) # Array containing the counters for permeation in each direction, according to NN (line 67 of the script).
permArray = numpy.hstack((numpy.linspace(0, len(traj), len(traj) + 1), permArray), dtype = int)

# def getOldPos():
#     atomIDs = pdb.select(selName).getIndices()
#     zetas = pdb[atomIDs].getCoords()[:,2] # Get z coordinate for all elements in the selection.
#     Arr = numpy.vstack((atomIDs, zetas))
#     return Arr

t1 = datetime.now()

if direction == 'up':
    counter = 0
    for i, frame in enumerate(traj):
        frame.superpose()
        posArray = atomsSelection.getCoords()
        flags = numpy.where([posArray[:,2] < 0, posArray[:,2] > 0], [numpy.ones(len(posArray)), 6*numpy.ones(len(posArray))], numpy.zeros(len(posArray)))
        # inCyl = numpy.argwhere((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < upZ) & (posArray[:,2] > loZ))
        inCyl = numpy.where(((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < upZ) & (posArray[:,2]) > loZ), 3*numpy.ones(len(posArray)), flags)

        diffArray = numpy.concatenate(numpy.argwhere(flags - flagsOld[:,-1] != 0))

        for value in diffArray:
            if 1 + numpy.sum((flagsOld[value], flags[value])) == 7:
                counter += 1
                permArray[i,-1] = counter
            numpy.roll(flagsOld[value],-1)
            flagsOld[value,-1] = flags[value]

        
elif direction == 'down':
    counter = 0
    for i, frame in enumerate(traj):
        frame.superpose()
        posArray = atomsSelection.getCoords()
        flags = numpy.where([posArray[:,2] < 0, posArray[:,2] > 0], [numpy.ones(len(posArray)), 6*numpy.ones(len(posArray))], numpy.zeros(len(posArray)))
        # inCyl = numpy.argwhere((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < upZ) & (posArray[:,2] > loZ))
        inCyl = numpy.where(((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < upZ) & (posArray[:,2]) > loZ), 3*numpy.ones(len(posArray)), flags)

        diffArray = numpy.concatenate(numpy.argwhere(flags - flagsOld[:,-1] != 0))

        for value in diffArray:
            if 1 + numpy.sum((flagsOld[value], flags[value])) == -3:
                counter += 1
                permArray[i, -1] = counter
            numpy.roll(flagsOld[value],-1)
            flagsOld[value,-1] = flags[value]

elif direction == 'both':
    counter = [0, 0]
    for i, frame in enumerate(traj):
        frame.superpose()
        posArray = atomsSelection.getCoords()
        flags = numpy.where([posArray[:,2] < 0, posArray[:,2] > 0], [numpy.ones(len(posArray)), 6*numpy.ones(len(posArray))], numpy.zeros(len(posArray)))
        # inCyl = numpy.argwhere((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < upZ) & (posArray[:,2] > loZ))
        inCyl = numpy.where(((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < upZ) & (posArray[:,2]) > loZ), 3*numpy.ones(len(posArray)), flags)

        diffArray = numpy.concatenate(numpy.argwhere(flags - flagsOld[:,-1] != 0))

        for value in diffArray:
            if 1 + numpy.sum((flagsOld[value], flags[value])) == 7:
                counter[0] += 1
            elif 1 + numpy.sum((flagsOld[value], flags[value])) == -3:
                counter[1] += 1

            permArray[-3] = counter[1] # Direction -z
            permArray[-2] = counter[0] # Direction +z
            permArray[-1] = counter[0] + counter[1] # Total count
            numpy.roll(flagsOld[value],-1)
            flagsOld[value,-1] = flags[value]

    else:
        print('PROBLEM: Direction incorrect!')
        exit()


permArray = permArray.astype('str')

outFile = open(outName, 'w+')
if direction == 'up':
    outFile.write('#Frame \t z+ permeations events \n')
elif direction == 'down':
    outFile.write('#Frame \t z- permeation events \n')
else:
    outFile.write('#Frame \t z- P. Events \t z+ P.Events \t Total \n')

for vals in permArray:
    outFile.write('{0} \n'.format(' \t'.join(vals)))
    
outFile.close()

t2 = datetime.now()

print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))

    # Ahora tengo que cambiar las flags de los índices que están dentro del cilindro SI ES QUE CORRESPONDE.

    # flags = numpy.where(z < loZ, numpy.ones(len(z)), numpy.zeros(len(z)))
    # flags = numpy.where(z > upZ, 4*numpy.ones(len(z)), flags)
    # flags = numpy.where(flags == 0, 2*numpy.ones(len(z)), flags)
    # for j in diffArray:
    #     for k in range(3):
    #         if dataDict[atomIDs[j]][k] != 0 and dataDict[atomIDs[j]][k] != flags[j]:
    #             dataDict[atomIDs[j]] = flags[j]



