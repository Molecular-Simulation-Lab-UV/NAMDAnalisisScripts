"""
This program relies on a user-defined cylinder that roughly defines the pore in order to count molecules
going from the inner solvent to the outer solvent through the pore. To achieve this, the variables
'upperZ', 'lowerZ' and 'rad' are needed in the input file, which differs from other scripts in this repository.
"""

import prody
import numpy
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = "Calculate the permeation events of selection's atoms through an equivalent cylinder.")
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

delta = 4
dcdName = []
outName = 'outFile.out'

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

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

# Setting up structure and linking the trajectory files to it.

pdb = prody.parsePDB(pdbName)
traj = prody.Trajectory(dcdName[0])
if len(traj) > 1:
    for d in dcdName[1:]:
        traj.addFile(d)

traj.link(pdb)
traj.setCoords(pdb)
traj.setAtoms(pdb.select(refName)) # refName = Selection used when aligning frames (frame.superpose())
atomsSelection = pdb.select(selName)
lenSel = len(atomsSelection)

NN = 3 # Counting in z+ direction, z- direction, and total events: that's the 3

atomIDs = atomsSelection.getIndices() # "sel" selection. Usually would be "name OH2".

"""
flags and flagsOld hold the indicators for signaling in which region the atom is.
1: z < 0 and outside the cylinder.
3: atom inside the cylinder.
6: z > 0 and outside the cylinder.
flags holds the indicator (flag) for each atom at the current frame/timestep.
flagsOld holds the previous 2 indicators - shape (len(traj), 2). Only changes if its [i-th, -1] value differs from the [i-th] value from the 'flags' array.
The corresponding element [frameNumber, -1] gets updated when it differs from the respective one from 'flags'.
permArray simply holds the frame and the counts in each direction, as well as the total.
diffArray holds the 'flags' ARRAY's index (not the atom's index!) of the differing flags.
"""
flagsOld = numpy.zeros((lenSel, 4))
permArray = numpy.zeros((len(traj), NN)) # Array containing the counters for permeation in each direction, according to NN (line 65 of the script)
permArray = numpy.hstack((numpy.linspace(1, len(traj), len(traj))[:,numpy.newaxis], permArray)).astype('int')

t1 = datetime.now()

counter = [0, 0]
for i, frame in enumerate(traj):
    prody.wrapAtoms(pdb, unitcell = frame.getUnitcell()[:3], center = prody.calcCenter(pdb.select(refName)))
    frame.superpose()
    posArray = pdb.select(selName).getCoords()
    # Flag atoms with postion z > 0 with 6, and z < 0 with 1
    flags = numpy.where(posArray[:,2] < 0, numpy.ones(lenSel), 5*numpy.ones(lenSel))
    # Include the flags for the atoms inside the cylinder
    flags = numpy.where(numpy.logical_and(posArray[:,2] < upZ - delta, posArray[:,2] > loZ + delta), 3*numpy.ones(lenSel), flags)
    flags = numpy.where((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < upZ) & (posArray[:,2] > upZ - delta), 4*numpy.ones(lenSel), flags)
    flags = numpy.where((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < loZ + delta) & (posArray[:,2] > loZ), 2*numpy.ones(lenSel), flags)
    mask = numpy.flatnonzero(numpy.logical_or(flags == 1, flags == 5))
    print(f'{flagsOld}   {flags}')
    try:
        diffArray = numpy.flatnonzero(numpy.all(numpy.not_equal(flagsOld[:, 1:], flags[:,numpy.newaxis]), axis = 1)) # Explained a while back.
    except:
        continue
    counter[1] += len(numpy.flatnonzero(numpy.all(numpy.equal(numpy.hstack((flagsOld[diffArray], flags[diffArray, numpy.newaxis])), [1,2,3,4,5]), axis = 1))) 
# Maybe directly modify the permArray instead, bypassing the counter variable.
    counter[0] += len(numpy.flatnonzero(numpy.all(numpy.equal(numpy.hstack((flagsOld[diffArray], flags[diffArray, numpy.newaxis])), [5,4,3,2,1]), axis = 1)))
        # Update the counters and flagsOld's flags.
    permArray[i,-3] = counter[0] # Direction -z
    permArray[i,-2] = counter[1] # Direction +z
    permArray[i,-1] = counter[0] + counter[1] # Total count
        # Push flagsOld flags (for the specific atomIndex) to make place for the updated one
        # flagsOld[value] = numpy.roll(flagsOld[value],-1)
        # flagsOld[value,-1] = flags[value] # Update flags in flagsOld
    flagsOld[diffArray] = numpy.roll(flagsOld[diffArray], -1)
    flagsOld[mask] = 0
    flagsOld[diffArray, -1] = flags[diffArray]
#    print(flagsOld[6792])
#    flagsOld[mask] = numpy.hstack((numpy.zeros((len(mask),3)), flags[mask, numpy.newaxis])) # Possible optimization here; the "try" statement could be shortened by avoiding the inclusion of the indices where flags == 5 or flags == 1.
#    print(f'Times: {dummyT2 - dummyT1}, {dummyT3 - dummyT2}')

permArray = permArray.astype('str') # Handy for writing to a file

outFile = open(outName, 'w+')
outFile.write('#Frame\t |\t z-\t |\t z+\t |\t Total Permeation Events\n')

for vals in permArray:
    outFile.write('{0} \n'.format(' \t '.join(vals)))
    
outFile.close()

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
