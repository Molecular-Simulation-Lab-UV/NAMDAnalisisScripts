#!/usr/bin/env python3

import prody
import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='Program to calculate distances between a first selection and all others')
parser.add_argument('-i', '--in_file', type=str, required=True, help='Path, either absolute or relative, to input file.')
parser.add_argument('-d', '--dimension', type=int, required=False, help='Distance measured: 1 = 1D, z-distance. 2 = 2D, sqrt(x^2 + y^2) distance. 3 = 3D distance. [Default]')
parser.add_argument('-c', '--cores', type=int, required=False, default=4, help='Core count assigned to the calculations. Default = 4')
parser.add_argument('-t', '--transform', required=False, action='store_true',
                     default=False, help='Whether or not to apply transformations to properly center the system. Default = False')

arg = parser.parse_args()
if arg.dimension == None:
    dim = 3
else:
    dim = arg.dimension

#  Should usually go lower, but we put it here to be able to select multiple dcds and selections
selName = {}
dcdName = []
mainSel = None

# Standard outFile name to avoid error when writing. If set in the input file, the name will be replaced.
outName = 'outFile.dat'

# =============================================================================
# Parsing the input file to gather inputs and selections.
# Ref is taken as an alignment selection. Then we search for sel. First in top-down arrangement is the main one.
# The distance is calculated between this first selection and all others.
# =============================================================================

print('\nReading input file and asigning variables')
inFile = open(arg.in_file, 'r')

for line in inFile:
    l = line.strip().split()
    if 'dcd' in l[0].lower():
        dcdName.append(l[1])
    elif l[0].lower() == 'ref':
        if len(l[1:]) > 1:
            refName = ' '.join(l[1:])
        else:
            refName = l[1]
    elif l[0].lower() == 'pdb':
        pdbName = l[1]
    elif l[0].lower() == 'psf':
        psfName = l[1]
    elif 'sel' in l[0].lower():
        if len(l[1:]) > 1 and mainSel == None:
            mainSel = ' '.join(l[1:])
        elif len(l[1:]) > 0 and mainSel != None:
            selName[l[0]] = ' '.join(l[1:])
        elif len(l[1:]) < 2 and mainSel == None:
            mainSel = l[1]
        else:
            selName[l[0]] = l[1]
    elif l[0].lower() == 'out':
        outName = l[1]

inFile.close()

universe = mda.Universe(psfName,
                        dcdName
                        )
reference = mda.Universe(psfName, pdbName)
alignment_ref = reference.select_atoms(refName)

ref = universe.select_atoms(refName)
atoms = universe.atoms

if arg.transform:
    workflow = [transformations.unwrap(atoms, max_threads=arg.cores),
                transformations.center_in_box(ref, center='geometry', max_threads=arg.cores),
                transformations.wrap(atoms, compound='segments', max_threads=arg.cores),
                transformations.fit_rot_trans(ref, alignment_ref, weights=None, max_threads=arg.cores)
    ]

    universe.trajectory.add_transformations(*workflow)


selections = [universe.select_atoms(mainSel)] + [universe.select_atoms(sel)
                                                  for sel in selName.values()]

t1 = datetime.now()
distArray = np.zeros((len(universe.trajectory), len(selName)))

def create_distance_calculator(dimension):
    if dimension == 1:
        def calcDist(centers):
            return centers[1:, -1] - centers[0, -1]
    
    elif dimension == 2:
        def calcDist(centers):
            diff = centers[1:, :2] - centers[0, :2]
            return np.linalg.norm(diff, axis=1)
    else:
        def calcDist(centers):
            diff = centers[1:] - centers[0]
            return np.linalg.norm(diff, axis=1)
    
    return calcDist

calcDist = create_distance_calculator(dim)

for f, frame in enumerate(universe.trajectory):
    centers = np.array([sel.center(None) for sel in selections])
    distArray[f] = calcDist(centers)

print('Done. Writing output file')

# For writing to outFile with "join" option.
distArray = distArray.astype('str')

outFile = open(outName, 'w+')
outFile.write('#Frame \t {0} \n'.format('\t '.join(selName.keys())))
for k, vals in enumerate(distArray):
    outFile.write('{0:7} \t {1} \n'.format(k, '\t'.join(vals)))
    
outFile.close()

t2 = datetime.now()
print('Done. Completion took {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))

###################################################
#################### LEGACY #######################
###################################################


# traj = prody.Trajectory(dcdName[0])
# pdb = prody.parsePDB(pdbName)

# # "Add" all of the trajectories, one after the other, together.
# if len(dcdName) > 1:
#     for i in dcdName[1:]:
#         traj.addFile(i)

# traj.link(pdb)
# traj.setCoords(pdb)
# traj.setAtoms(pdb.select(refName))

# # Calculate Z distance between selections
# # This function returns the absolute value. You'll have to change the "return" at
# # the end of the following function in order to get the vector component.
# def calc1D(sel1, sel2):
#     selPos1 = prody.calcCenter(sel1)[2]
#     selPos2 = prody.calcCenter(sel2)[2]
#     return abs(selPos1 - selPos2)
    
# # Calculate radial distance between selections
# def calc2D(sel1, sel2):
#     selPos1 = prody.calcCenter(sel1)[:2]
#     selPos2 = prody.calcCenter(sel2)[:2]
#     return numpy.linalg.norm(selPos1 - selPos2)

# # Set the first element as the one which all distance measurements will be
# # calculated against.

# refSel = pdb.select(mainSel)
# print(refSel)

# selections = []
# for key in selName:
#     selections.append(pdb.select(selName[key]))

# print('\nBeginning distance calculations for {0} frames'.format(len(traj)))
# t1 = datetime.now()

# # Calculate the distances
# distArray = numpy.zeros((len(traj), len(selName)))
# print(distArray.shape)

# if dim == 1:
#     for i, frame in enumerate(traj):
#         frame.superpose()
#         for j, sel in enumerate(selections):
#             distArray[i,j] = calc1D(refSel, sel)
# elif dim == 2:
#     for i, frame in enumerate(traj):
#         frame.superpose()
#         for j, sel in enumerate(selections):
#             distArray[i,j] = calc2D(refSel, sel)
# elif dim == 3:
#     for i, frame in enumerate(traj):
#         frame.superpose()
#         for j, sel in enumerate(selections):
#             distArray[i,j] = prody.calcDistance(prody.calcCenter(refSel), prody.calcCenter(sel))

# else:
#     print('Something went wrong. Check dimension (-d) required (╯°□°)╯︵ ┻━┻')

# print('Done. Writing output file')

# # For writing to outFile with "join" option.
# distArray = distArray.astype('str')

# outFile = open(outName, 'w+')
# outFile.write('#Frame \t {0} \n'.format('\t '.join(selName.keys())))
# for k, vals in enumerate(distArray):
#     outFile.write('{0:7} \t {1} \n'.format(k, '\t'.join(vals)))
    
# outFile.close()

# t2 = datetime.now()
# print('Done. Completion took {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))