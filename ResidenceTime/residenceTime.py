import prody
import argparse
from datetime import datetime
import numpy

parser = argparse.ArgumentParser(description = 'Calculate the collective dipole of the given selection.')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

dcdName = []
outName = 'outFile.out'

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
    if len(l) > 0:
        if 'dcd' in l[0]: # Any amount of .dcd files
            dcdName.append(l[1])
        elif l[0] == 'pdb':
            pdbName = l[1]
        elif l[0] == 'ref': # Alignment selection
            if len(l[1:]) > 1:
                refName = ' '.join(l[1:])
            else:
                refName = l[1]
        elif 'sel' in l[0]: # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2"
            if len(l[1:]) > 1:
                selName = ' '.join(l[1:])
            else:
                selName = l[1]
        elif l[0] == 'out':
            outName = l[1]
        if l[0].lower() == 'zmin':
            zMin = float(l[1])
        elif l[0].lower() == 'zmax':
            zMax = float(l[1])
        elif l[0].lower() == 'nbins':
            nBins = int(l[1])
        elif l[0].lower() == 'rad':
            rad = float(l[1])
        elif l[0].lower() == 'thr':
            thr = int(l[1]) # Can be a space-separated list. Each element is a number of minimum frames that define a "time-wise bin"
                            # E.g. If 5, 100, 500 is given, data about residence time in intervals < 5 | 5 - 100 | 100 - 500 | 500 > will be returned.
                            # The first interval is not considered. If set to 0, all actual intervals are considered.
    else:
        pass

inFile.close()

pdb = prody.parsePDB(pdbName)
dcd = prody.Trajectory(dcdName[0])
if len(dcdName) > 1:
    for d in dcdName[1:]:
        dcd.addFile(d)

dcd.link(pdb)
dcd.setCoords(pdb)
dcd.setAtoms(pdb.select(refName)) # refName = Selection used when aligning frames (frame.superpose())
atomsSelection = pdb.select(selName)

# The number of frames that a molecule has stayed in one same bin
lCyl = abs(zMax) + abs(zMin)
binSize = lCyl/nBins
binArray = numpy.arange(zMin, zMax + binSize, binSize)

t1 = datetime.now()

moleculesInBins = numpy.zeros((len(atomsSelection), nBins+1)).astype(int)
moleculesInBins[:,0] = atomsSelection.getIndices()
oldMoleculesInBins = moleculesInBins.copy()
# The value for each bin throughout the simulation.
binsInTime = numpy.zeros((len(dcd), nBins))
f0 = dcd.next()
prody.wrapAtoms(pdb, unitcell = f0.getUnitcell()[:3], center = prody.calcCenter(pdb.select(refName)))
f0.superpose()
oldPos = pdb.select(f'{selName} and (x^2 + y^2) < {rad**2} and z > {zMin} and z < {zMax}').getCoords()[:,-1]
oldInd = pdb.select(f'{selName} and (x^2 + y^2) < {rad**2} and z > {zMin} and z < {zMax}').getIndices()
oldInBin = numpy.argwhere((oldPos[:,numpy.newaxis] >= binArray[numpy.newaxis,:-1]) & (oldPos[:,numpy.newaxis] < binArray[numpy.newaxis,1:]))
# Initializers for positions, indices and whether it's in bin or not at frame 0
for f, frame in enumerate(dcd):
        prody.wrapAtoms(pdb, unitcell = frame.getUnitcell()[:3], center = prody.calcCenter(pdb.select(refName)))
        frame.superpose()
        sel = pdb.select(f'{selName} and (x^2 + y^2) < {rad**2} and z > {zMin} and z < {zMax}')
        pos = sel.getCoords()[:,-1] # Grab selection position, z-coordinate
        ind = sel.getIndices() # Grab selection (atom) indices
        # Defines which INDEX of the pos/ind array goes into which bin.
        inBin = numpy.argwhere((pos[:,numpy.newaxis] >= binArray[numpy.newaxis,:-1]) & (pos[:,numpy.newaxis] < binArray[numpy.newaxis,1:]))
        # Combine the index of an atom in the selection with the respective bin.
        horizontalStack = numpy.vstack((ind, inBin[:,1] + 1)).T
        oldHorizontalStack = numpy.vstack((oldInd, oldInBin[:,1] + 1)).T
        # Compare which pair index/bin persists after the previous frame
        mask = numpy.flatnonzero((oldHorizontalStack == horizontalStack[:,None]).all(-1).any(-1))
        # Check which of the global selection atoms have remained in their bin (moleculesMask) and which haven't (notMoleculesMaks)
        moleculesMask = numpy.isin(moleculesInBins[:,0], ind[mask])
        notMoleculesMask = numpy.invert(moleculesMask)
        # Add +1 to the counter for those atoms that persist
        moleculesInBins[moleculesMask, inBin[mask,1]+1] += 1
        # Select those atoms that 1. Haven't remained in their bin and 2. Are higher than the threshold established.
        validTimes = numpy.where((moleculesInBins[notMoleculesMask, 1:] > thr), moleculesInBins[notMoleculesMask, 1:], 0).astype(float)
        nonZeroBins = numpy.count_nonzero(validTimes, axis=0)
        # Add the average counter, of all atoms, per bin for the current frame
        binsInTime[f] = numpy.where(nonZeroBins != 0, validTimes.sum(axis=0)/nonZeroBins, 0)
        # Reset the counter for those that don't persist
        moleculesInBins[notMoleculesMask, 1:] = 0
        # Update values
        oldPos = pos
        oldInd = ind
        oldInBin = inBin

# Get the average in the time axis and output
totalFrames = numpy.sum(binsInTime, axis = 0)
nonZeroFrames = numpy.count_nonzero(binsInTime, axis = 0)
finalArray = numpy.where(nonZeroFrames != 0, totalFrames / nonZeroFrames, 0)
outFile = open(outName, 'w+')
outFile.write('# Z |\t Residence time avg \n')

for f, vals in enumerate(finalArray):
    outFile.write('{0} \t {1} \n'.format(binArray[f+1] - binSize/2, vals))
    
outFile.close()

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))
