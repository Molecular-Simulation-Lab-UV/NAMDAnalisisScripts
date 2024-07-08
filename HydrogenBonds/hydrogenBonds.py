from symbol import atom
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis import transformations
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = 'Calculate the collective dipole of the given selection.')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')
parser.add_argument('-u', '--update', action = 'store_true', required = False, default = False, help = 'Calculate the average for the trajectory if included (with "-a" option), frame-by-frame-wise if not included. Default: False')
parser.add_argument('-b', '--bins', action = 'store_true', required = False, default = False, help = 'Return dipole of selection "sel" binned along a cylinder. Default: False')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

dcdName = []
outName = 'outFile.out'

for line in inFile:
    l = line.strip().split()
    if len(l) > 0:
        if 'dcd' in l[0]: # Any amount of .dcd files
            dcdName.append(l[1])
        elif l[0] == 'pdb':
            pdbName = l[1]
        elif l[0] == 'psf':
            psfName = l[1]
        elif l[0] == 'ref': # Alignment selection
            if len(l[1:]) > 1:
                refName = ' '.join(l[1:])
            else:
                refName = l[1]
        elif 'se1' in l[0]: # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2"
            if len(l[1:]) > 1:
                selName1 = ' '.join(l[1:])
            else:
                selName1 = l[1]
        elif 'sel2' in l[0]: # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2"
            if len(l[1:]) > 1:
                selName2 = ' '.join(l[1:])
            else:
                selName2 = l[1]
        elif l[0] == 'out':
            outName = l[1]
        if arg.bins == True:
            if l[0].lower() == 'zmin':
                zMin = float(l[1])
            elif l[0].lower() == 'zmax':
                zMax = float(l[1])
            elif l[0].lower() == 'nbins':
                nBins = int(l[1])
            elif l[0].lower() == 'rad':
                rad = float(l[1])
        elif arg.bins == False and warning:
            print("\nWARNING: --bins (-b) flag was not set. Data won't be binned.\n")
            warning = False
    else:
        pass

inFile.close()

# refU = mda.Universe(psfName, pdbName)
u = mda.Universe(psfName, dcdName)
# sel = u.select_atoms(refName)
# refSel = refU.select_atoms(refName)
# atomGroup = u.atoms

# workflow = [transformations.unwrap(atomGroup),
#             transformations.center_in_box(sel, center = 'mass'),
#             transformations.unwrap(atomGroup, compound = 'fragments'),
#             transformations.fit_rot_trans(sel, refSel)]
# u.trajectory.add_transformations(*workflow)

hbonds = HBA(universe = u, hydrogens_sel = selName1, acceptors_sel = selName2, d_a_cutoff = 3.0, d_h_a_angle_cutoff = 140, update_selections = arg.update)
hbonds.run(verbose = True)

