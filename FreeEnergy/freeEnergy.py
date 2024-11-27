# 
import linecache
from collections import namedtuple
import re

import numpy as np
import matplotlib.pyplot as plt
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='Free Energy calculation from energy state frequency')
parser.add_argument('-i', '--in_file', type=str, required=True, help='Path, either absolute or relative, to input file.')

arg = parser.parse_args()

colvars = {}
colvar_config = {}
dcdName = []
mainSel = None
R = 8.314
T = 298.15

inFile = open(arg.in_file, 'r')

def read_namd_cv(cv, cv_dict):
    headers = list(map(str.strip, linecache.getline(cv[1], 1).split()[1:]))
    data = np.loadtxt(cv[1])
    idx = headers.index(cv[0])
    cv_dict[headers[idx] + '_namd'] = data[:,idx]

def parse_colvar_bounds(filepath):
    """
    Parse a configuration file to extract colvar parameters.
    Returns a list of named tuples containing lowerBoundary, upperBoundary, and width
    for each colvar block found.
    
    Args:
        filepath (str): Path to the configuration file
        
    Returns:
        list[ColvarBounds]: List of named tuples containing the extracted values
    """
    # Define the named tuple structure
    Colvar = namedtuple('colvar', ['name', 'lower_boundary', 'upper_boundary', 'width'])
    
    # Initialize results list
    results = []
    
    # Read the file content
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Find all colvar blocks using regex
    colvar_pattern = r'(?i)colvar\s*{([^}]+)}'
    colvar_blocks = re.finditer(colvar_pattern, content)
    
    for block in colvar_blocks:
        block_content = block.group(1)
        
        # Extract name
        name_match = re.search(r'name\s+(\w+)', block_content)
        name = name_match.group(1) if name_match else None
        
        # Extract values using regex
        lower = re.search(r'lowerBoundary\s+([-\d.]+)', block_content)
        upper = re.search(r'upperBoundary\s+([-\d.]+)', block_content)
        width = re.search(r'width\s+([-\d.]+)', block_content)
        
        # Convert to float if found, None if not found
        # TODO: Replace with true default values
        lower_val = float(lower.group(1)) if lower else None
        upper_val = float(upper.group(1)) if upper else None
        width_val = float(width.group(1)) if width else None
        
        # Create named tuple and add to results
        if any([lower_val, upper_val, width_val]):  # Only add if at least one value was found
            bounds = Colvar(
                name=name,
                lower_boundary=lower_val,
                upper_boundary=upper_val,
                width=width_val
            )
            results.append(bounds)
    
    return results

for line in inFile:
    l = line.strip().split()
    if l[0].lower() == 'cv_file':
        cv_file = l[1]
    elif 'cv_namd' in l[0].lower():
        cv = l[1].split(':')
        read_namd_cv(cv, colvars)
    elif 'cv' in l[0].lower():
        cv = l[1].split(':')
        dummy_dat = np.loadtxt(cv[1])
        colvars[cv[0]] = dummy_dat[:,-1]
    elif l[0] == 'pdb':
        pdbName = l[1]
    elif l[0] == 'out':
        outName = l[1]

inFile.close()

colvar_config = parse_colvar_bounds(cv_file)
fes_grid = [np.arange(cv_i.lower_boundary, cv_i.upper_boundary + cv_i.width, cv_i.width) for cv_i in colvar_config]
stacked_colvars = np.vstack(list(colvars.values())).T
# Create grid and bins for arbitrary dimensions
# stacked_colvars should now be [60001, N] where N is the number of colvars
counts, _ = np.histogramdd(
    stacked_colvars,  # Data points
    bins=fes_grid     # Bins for each dimension
)
# np.savetxt('./test_counts.txt', counts, fmt='%.2e')
# TODO: Generate figure only if 1D or 2D
counts[counts==0] = 1
fes = -R * T * np.log(counts/stacked_colvars.shape[0])
X, Y = np.meshgrid(*fes_grid)
X = X[:-1, :-1] + (X[1:, 1:] - X[:-1, :-1])/2
Y = Y[:-1, :-1] + (Y[1:, 1:] - Y[:-1, :-1])/2
plt.figure(figsize=(10,7))
plt.contourf(X.T, Y.T, fes, 20, cmap='rainbow')
plt.colorbar()
plt.savefig('test.png', format='png')
