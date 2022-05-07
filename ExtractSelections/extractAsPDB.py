"""
Batch isolate 'selection' from all the files found in -p path. These files have to be .pdbs.
Outputs have the same name, with at most the 2 first words in 'selection' appended to the file name.
The main purpose is to simplify the process of generating elements en masse for, let's say, docking tests.
"""

import prody
import argparse
import glob

# Default selection.
sel = "protein"

parser = argparse.ArgumentParser(description = 'Extracts only the selection passed as argument -s from a list of pdbs found in the input path. Then writes them as another .pdb with the same name, plus the selection name (up to 2 words to avoid cluttering), to the same directory.')
parser.add_argument('-p', '--path', type = str, required = True, help = 'Path, either absolute or relative, to folder containing pdbs')
parser.add_argument('-s', '--selection', type = str, required = False, help = '[Optional] VMD-like selection to be isolated and written as an independent .pdb. Must use "" when passing the argument (e.g.: "resid 20"). Defaults to "protein".')

arg = parser.parse_args()
if arg.selection != None:
	sel = arg.selection
fileDir = glob.glob(arg.path + '/*.pdb')
outputs = []
for file1 in fileDir:
	f1 = file1.strip().split('.pdb')[0]
	outputs.append(f1)

for i in range(len(fileDir)):
	pdb = prody.parsePDB(fileDir[i])
	molecule = pdb.select(sel).copy()
	if len(sel.split()) >= 2:
		prody.writePDB(outputs[i] + '.' + ''.join(sel.split()[:2]) + '.pdb', molecule)
	else:
		prody.writePDB(outputs[i] + '.' + sel + '.pdb', molecule)
print('{} molecules extracted from path {}'.format(len(fileDir), arg.path))
