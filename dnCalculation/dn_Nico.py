#################
#### IMPORTS ####
#################

from datetime import datetime
import prody
import numpy as np
import argparse
import dnModule_Nico as dn

###################
#### VARIABLES ####
###################

parser = argparse.ArgumentParser(description = "Calculate the permeation events of selection's atoms through a pore.")
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
    if 'dcd' in l[0]: # Any amount of .dcd files. Multiple lines, in order up-down, can be input. Only requisite is "dcd" should appear as keyword.
        dcdName.append(l[1])
    elif l[0].lower() == 'pdb':
        pdbName = l[1]
    elif l[0].lower() == 'ref': # Alignment selection for superposition of the .dcd frames onto the pdb conformation.
        if len(l[1:]) > 1:
            refName = ' '.join(l[1:])
        else:
            refName = l[1]
    elif l[0].lower() == 'rad': # Radius of the cylinders operating as checkpoints in the entry and exit of the pore.
        rad = float(l[1])
        rad2 = rad**2
    elif l[0].lower() == 'nbins': # Number of bins for cylinder subdivision.
        binNumber = float(l[1])
    elif l[0].lower() == 'upperz': # Upper z coordinate/position of the cylinder.
        upZ = float(l[1])
    elif l[0].lower() == 'lowerz': # Lower z coordinate/position of the "pore".
        loZ = float(l[1])
    elif l[0].lower() == 'sel': # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if len(l[1:]) > 1:
            selName = ' '.join(l[1:])
        else:
            selName = l[1]
    elif l[0].lower() == 'out': # Path to the output file, name included.
        outName = l[1]

lCyl = abs(upZ) + abs(loZ)
binSize = lCyl/binNumber

###################
#### ARCHIVOS  ####
###################

print("zmax: ", upZ)
print("zmin: ", loZ)
print("radio: ", rad)
print("largo del cilindro: ", lCyl)
print("cantidad de bins: ", binNumber)
print("tamaño de bins: ", binSize)

t1 = datetime.now()

traj = prody.Trajectory(dcdName[0])
for d in dcdName[1:]:
    traj.addFile(d)

pdb = prody.parsePDB(pdbName) #Cargamos la información de coordenadas del pdb.

traj.setAtoms(pdb.select(refName)) #Alineamos respecto al pdb de referencia usando la selección especificada.
traj.setCoords(pdb) #Establecemos el pdb como frame de referencia
traj.link(pdb)

dn = dn.dnMatrixCalculation(traj, pdb, loZ, upZ, rad2, binSize, refName)
np.save(outName, dn)

t2 = datetime.now()

print('\nFIN')
print('Time to completion was {0}'.format((t2.replace(microsecond=0) - t1.replace(microsecond=0))))