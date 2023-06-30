"""
This is mainly for my own use and it was created because I'm lazy c:
It's not warranted to work... at all. Use at your own risk.
"""


import numpy
import argparse
from glob import glob
from os import getcwd

parser = argparse.ArgumentParser(description = "Concatenate two or more simple text files; made for time series data.")
parser.add_argument('-o', '--outfile', type = str, required = False, help = 'Output file path and name. Default is output.out.')
parser.add_argument('-i', '--in_file', type = str, required = True, nargs = '+', help = 'List of files to join. Wildcards (*) can be used.')
parser.add_argument('-c', '--columns', type = int, default = 0, nargs = '+', required = False, help = "Column indices to modify so they start from previous file's last step.\nDefault is only first (column 0).")

outName = 'output.out' # Default output file name

args = parser.parse_args()

if '*' not in args.in_file:
    files = args.in_file
else:
    print('* wildcard is not adviceable, as the built-in sorted function will be used to order the files.')
    files = sorted(glob(getcwd() + '/' + args.in_file))

data = []
if args.columns is not None:
    columns = args.columns
else:
    columns = [0]

if args.outfile is not None:
    outName = args.outfile

data = []

# First file, not to be modified
with open(files[0], 'r') as f1:
    for line in f1:
        firstLine = line
        break
dat = numpy.loadtxt(files[0])
data.append(dat)
lastLine = dat[-1, columns]

# Remaining files
for file in files[1:]:
    dat  = numpy.loadtxt(file)
    dat[:, columns] = dat[:, columns] + lastLine
    lastLine = dat[-1, columns]
    data.append(dat)

data = numpy.vstack(data)
data = data.astype('str') # Handy for writing to a file

outFile = open(outName, 'w+')
outFile.write(firstLine)

for vals in data:
    outFile.write('{0} \n'.format(' \t '.join(vals)))
    
outFile.close()

print('\nFIN')