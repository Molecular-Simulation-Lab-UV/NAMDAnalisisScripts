import prody
import numpy
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description = 'Calculation of the dn variable, according to Zhu and Zon 2004')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, either absolute or relative, to the input file.')