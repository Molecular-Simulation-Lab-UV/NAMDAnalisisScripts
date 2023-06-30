import numpy
import argparse
from glob import glob
from os import getcwd

parser = argparse.ArgumentParser(description = "Concatenate two or more simple text files; made for time series data.")
parser.add_argument('-i', '--in_file', type = str, required = True, nargs = '+', help = 'List of files to join. Wildcards (*) can be used.')
parser.add_argument('-s', '--simple', type = int, required = True, help = 'Simply join the files, or make one continue where the previous finished (for time series data).')
parser.add_argument('-c', '--columns', type = int, default = 0, nargs = '+', required = False, help = "Column indices to modify so they start from previous file's last step.\nDefault is only first (column 0).")
parser.add_argument('-o', '--outfile', type = int, required = False, help = 'Output file path and name. Default is output.out.')

outName = 'output.out' # Default output file name

args = parser.parse_args()

if '*' not in args.in_file:
    files = args.in_file
else:
    print('* wildcard is not adviceable, as the built-in sorted function will be used to order the files.')
    files = sorted(glob(getcwd() + '/' + args.in_file))

data = []
for file in files:
    try:
        dat  = numpy.loadtxt(file)
        lastLine = dat[-1]
    except:
        f1 = open(file, 'r')
