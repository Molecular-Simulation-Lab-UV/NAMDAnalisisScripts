"""
This program sets up two cylinders, a lower one and an upper one. The cylinders go from
'lowerZ' to 'lowerZ' + delta, and from 'upperZ' - delta to 'upperZ'; delta defaults to the value 3 when it's not specified.
Each cylinder represents an entry or exit of the pore when the protein's pore is aligned with the z axis. 
The radius of each cylinder is 'rad', as read from the input file. The script separates the simulation volume in 5 zones:
- A lower z zone:       for all atoms below z = lowerZ
- An upper z zone:      for all atoms above z = upperZ
- An intermediate zone: for all atoms in the range lowerZ + delta < z < upperZ - delta
- A lower cylinder:     for all atoms inside the cylinder of radius 'rad' and lowerZ < z < lowerZ + delta
- An upper cylinder:    for all atoms inside the cylinder of radius 'rad' and upperZ - delta < z < upperZ
Therefore, selection of the delta value should consider timestep of the .dcd file, which could make it so
molecules/atoms passing through the pore are not seen inside the cylinder due to a too large time-skip; and
also pore size, both length and approximate radius. Since the cylinder approximation is only implemented in the entry/exit points,
One can adjust the radius by measuring the distance of selected atoms in the entry/exit point.
"""

import prody
import numpy
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = "Calculate the permeation events of selection's atoms through a pore.")
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

delta = 3
dcdName = []
outName = 'outFile.out'

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
    if 'dcd' in l[0]: # Any amount of .dcd files. Multiple lines, in order up-down, can be input. Only requisite is "dcd" should appear as keyword.
        dcdName.append(l[1])
    elif l[0] == 'pdb':
        pdbName = l[1]
    elif l[0] == 'ref': # Alignment selection for superposition of the .dcd frames onto the pdb conformation.
        if len(l[1:]) > 1:
            refName = ' '.join(l[1:])
        else:
            refName = l[1]
    elif l[0].lower() == 'rad': # Radius of the cylinders operating as checkpoints in the entry and exit of the pore.
        rad = float(l[1])
        rad2 = rad**2
    elif l[0].lower() == 'delta': # Height of the cylinders operating as checkpoints in the entry and exit points of the pore. Defaults to 3 if not specified.
        delta = float(l[1])
    elif l[0].lower() == 'upperz': # Upper z position of the "pore" itself, as is considered by the user.
        upZ = float(l[1])
    elif l[0].lower() == 'lowerz': # Lower z position of the "pore".
        loZ = float(l[1])
    elif 'sel' in l[0]: # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if len(l[1:]) > 1:
            selName = ' '.join(l[1:])
        else:
            selName = l[1]
    elif l[0] == 'out': # Path to the output file, name included.
        outName = l[1]

inFile.close()

# Setting up structure and linking the trajectory files to it.

pdb = prody.parsePDB(pdbName)
traj = prody.Trajectory(dcdName[0])
if len(dcdName) > 1:
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
flags and flagsOld hold the indicators for signaling in which region each atom in the selection is.
1: z < 0 and outside the pore.
2: lowerZ < z < lowerZ + delta. Inside the lower cylinder, which acts as entry/exit checkpoint.
3: atom inside the pore, but not inside any of the cylinders. This assumes there's no radial way to exit the pore.
4: upperZ - delta < z < upperZ. Inside the upper cylinder, which also acts as exit/entry checkpoint.
5: z > 0 and outside the pore.
flags holds the indicator (flag) for each atom at the current frame/timestep.
flagsOld holds the previous 4 indicators - shape (len(traj), 4). Only changes if it detects a new value in flags that is not contained in flagsOld, at that specific array's index.
If flags[i] is either 1 or 5, the array resets that line to [0,0,0,1] or [0,0,0,5], respectively.
However, if flags[i] is 1 or 5 and flagsOld[i, :] is [5,4,3,2] or [1,2,3,4], respectively, we consider that a permeation event and increase the corresponding counter by 1.
permArray simply holds the frame number and the counts in each direction, as well as the total, at each frame.
diffArray holds the 'flags' ARRAY's index (not the atom's index!) of the differing flags; i.e, those indices that are not present in flagsOld.
"""
flagsOld = numpy.zeros((lenSel, 4)) # Initializing a zeros array for flagsOld
permArray = numpy.zeros((len(traj), NN)) # Array containing the counters for permeation in each direction (up, down, total)
permArray = numpy.hstack((numpy.linspace(1, len(traj), len(traj))[:,numpy.newaxis], permArray)).astype('int') # Pre-filled with the frame number.

t1 = datetime.now()

counter = [0, 0] # [up, down] permeation events counters
for i, frame in enumerate(traj):
    prody.wrapAtoms(pdb, unitcell = frame.getUnitcell()[:3], center = prody.calcCenter(pdb.select(refName)))
    frame.superpose()
    # Get selection positions.
    posArray = pdb.select(selName).getCoords()
    # Flag atoms with position z > 0 with 5, and z < 0 with 1. Those are the atoms outside the pore
    flags = numpy.where(posArray[:,2] < 0, numpy.ones(lenSel), 5*numpy.ones(lenSel))
    # Include the flags for the atoms inside the cylinders and between them. It's like a sandwich; a slab in between two cylinders
    flags = numpy.where(numpy.logical_and(posArray[:,2] < upZ - delta, posArray[:,2] > loZ + delta), 3*numpy.ones(lenSel), flags) # Slab
    flags = numpy.where((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < upZ) & (posArray[:,2] > upZ - delta), 4*numpy.ones(lenSel), flags) # Upper cylinder
    flags = numpy.where((numpy.sqrt(posArray[:,0]**2 + posArray[:,1]**2) < rad) & (posArray[:,2] < loZ + delta) & (posArray[:,2] > loZ), 2*numpy.ones(lenSel), flags) # Lower cylinder
    mask = numpy.flatnonzero(numpy.logical_or(flags == 1, flags == 5)) # Getting the indices of the flags array's elements that are outside the pore.
#    print(f'{flagsOld}   {flags}') A remnant from desperate times.

    try:
        diffArray = numpy.flatnonzero(numpy.all(numpy.not_equal(flagsOld[:,1:], flags[:,numpy.newaxis]), axis = 1)) # Explained a while back. Grabbing the differing indices to update them in flagsOld
    except:
        continue
    counter[1] += len(numpy.flatnonzero(numpy.all(numpy.equal(numpy.hstack((flagsOld[diffArray], flags[diffArray, numpy.newaxis])), [1,2,3,4,5]), axis = 1))) # Update up direction counter
# Maybe directly modify the permArray instead, bypassing the counter variable?
    counter[0] += len(numpy.flatnonzero(numpy.all(numpy.equal(numpy.hstack((flagsOld[diffArray], flags[diffArray, numpy.newaxis])), [5,4,3,2,1]), axis = 1))) # Update down direction counter
# Update the final array, the one we will output.
    permArray[i,-3] = counter[0] # Direction -z
    permArray[i,-2] = counter[1] # Direction +z
    permArray[i,-1] = counter[0] + counter[1] # Total count
# Push flagsOld flags (for the specific index from diffArray) to make place for the updated one
    flagsOld[diffArray] = numpy.roll(flagsOld[diffArray], -1)
# If the atom is outside the pore, reset its flags to 0
    flagsOld[mask] = 0
# Update all the latest flags in flagsOld. This allows us to check a history of the zones the atom has visited.
# This history is reset/deleted if the atom is found outside the pore, but it already was counted as a permeation event if applicable.
    flagsOld[numpy.union1d(diffArray, mask), -1] = flags[numpy.union1d(diffArray, mask)]

permArray = permArray.astype('str') # Handy for writing to a file

outFile = open(outName, 'w+')
outFile.write('#Frame\t |\t z-\t |\t z+\t |\t Total Permeation Events\n')

for vals in permArray:
    outFile.write('{0} \n'.format(' \t '.join(vals)))
    
outFile.close()

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
