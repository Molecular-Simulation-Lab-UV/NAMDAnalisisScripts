import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = "Calculate the permeation events of selection's atoms through a pore.")
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')

dcdName = []
mainSel = None
selName = None
outName = 'outFile.out'

# Reading inputs from the input parameters file. This is for the analysis calculation,
# NOT for the simulation itself.

print('\nReading the input and asigning variables')

for line in inFile:
    l = line.strip().split()
    if 'dcd' in l[0].lower(): # Any amount of .dcd files. Multiple lines, in order up-down, can be input. Only requisite is "dcd" should appear as keyword.
        dcdName.append(l[1])
    elif l[0].lower() == 'pdb':
        pdbName = l[1]
    elif l[0].lower() == 'ref': # Alignment selection for superposition of the .dcd frames onto the pdb conformation.
        if len(l[1:]) > 1:
            refName = ' '.join(l[1:])
        else:
            refName = l[1]
    elif 'sel' in l[0].lower(): # Should be a broader selection, not including the position bounds (without "and z <= XX") E.g., "name OH2" or "index 9000".
        if mainSel == None:
            if len(l[1:]) > 1:
                mainSel = ' '.join(l[1:])
            else:
                mainSel = l[1]
        elif mainSel != None and selName == None:
            if len(l[1:]) > 1:
                selName = ' '.join(l[1:])
            else:
                selName = l[1]
        else:
            print('Only two selections are needed to calculate RDF.')
            print('Selection {0} will be ignored.'.format(' '.join(l[1:])))
            continue
    elif l[0].lower() == 'out': # Path to the output file, name included.
        outName = l[1]

inFile.close()

universe = mda.Universe(pdbName, dcdName)
sel1 = universe.select_atoms(mainSel)
sel2 = universe.select_atoms(selName)

rdf_run = InterRDF(sel1, sel2)
rdf_run.run()

res = rdf_run.results
bins = res.bins
rdf = res.rdf
counts = res.count.astype(int)

outFile = open(outName, 'w+')
outFile.write('# Bin \t |\t RDF \t |\t Count \n')

for f in range(len(rdf)):
    outFile.write('{0:5.2f} \t {1:10.7f} \t {2:7d} \n'.format(bins[f], rdf[f], counts[f]))
    
outFile.close()
