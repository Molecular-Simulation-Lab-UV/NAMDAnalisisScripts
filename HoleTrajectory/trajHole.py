"""
Script for analyzing a dcd trajectory with the HOLE program. Said program (HOLE) has to be downloaded separately.
The output of this script is a file with averages, per bin, of the radius of the pore, and the asociated standard deviations.
First column: z-position of the middle of the bin.
Second column: average radius corresponding to that bin.
Third column: standard deviation for that bin.
"""

import numpy
import argparse
import MDAnalysis
import MDAnalysis.analysis.hole2 as hole2

holePath = '/opt/hole/hole2/exe/hole' # It's supposed to detect the binary in the path (module load hole2/hole), but never worked for me.
outName = 'outFile.dat'
binParams = [100, -20, 20]
delFiles = True

parser = argparse.ArgumentParser(description = 'Run HOLE and prepare an output file ready for plotting.\n File should contain at least a "dcd" and "psf" definition, and "bin".')
parser.add_argument('-i', '--in_file', type = str, required = True, help = 'Path, absolute or relative, to input file.')
parser.add_argument('-p', '--hole_path', type = str, required = False, help = 'Path to the HOLE binary. Can be set explicitly in the script as variable (check line 14 of this script.')

arg = parser.parse_args()
inFile = open(arg.in_file, 'r')
for line in inFile:
	l = line.strip().split()
	if l[0] == 'dcd':
		dcdName = l[1]
	elif l[0] == 'psf':
		psfName = l[1]
	elif dcdName != None and psfName != None and l[0] == 'sel':
		selName = ' '.join(l[1:])
	elif l[0] == 'out':
		outName = l[1]
	elif l[0] == 'keep':
		delFiles = False
	elif l[0] == 'bin':
		binParams = [int(l[1]), float(l[2]), float(l[3])] # l[1] should be number of bins, l[2] lower bound of pore Length, l[3] upper bound.
	
inFile.close()

if arg.hole_path != None:
	holePath = arg.hole_path

uni = MDAnalysis.Universe(psfName, dcdName)
if 'selName' in locals() and len(selName) >= 1:
	traj = hole2.HoleAnalysis(uni, select = selName, executable = holePath)

else:
	traj = hole2.HoleAnalysis(uni, executable = holePath)

traj.run(random_seed = 1000)


# Data ready for plot, adjusted to the stolen plot code

try:
	binned, bins = traj.bin_radii(bins = binParams[0], range = (binParams[1], binParams[2]))
except:
	print('Try and reduce the cylinder upper and lower z-values')

mean = numpy.array(list(map(numpy.mean, binned)))
midpoints = 0.5*bins[1:] + bins[:-1]
std = numpy.array(list(map(numpy.std, binned)))


outFile = open(outName, 'w+')
outFile.write('#Z position \t Mean radius value \t Standard deviation \n')
for i in range(len(mean)):
	outFile.write('{0:10.3f} \t {1:10.3f} \t {2:10.3f} \n'.format(midpoints[i], mean[i], std[i]))

outFile.close()
if delFiles == True:
	traj.delete_temporary_files()


#	Robado de HoleAnalysis (MDAnalysis --> analysis --> hole.py)
#        binned, bins = self.bin_radii(frames=frames, bins=bins, range=range)
#        mean = np.array(list(map(np.mean, binned)))
#        midpoints = 0.5 * bins[1:] + bins[:-1]
#
#        fig, ax = plt.subplots()
#        if n_std:
#            std = np.array(list(map(np.std, binned)))
#            ax.fill_between(midpoints, mean-(n_std*std), mean+(n_std*std),
#                            color=color, alpha=fill_alpha,
#                            label='{} std'.format(n_std))
#        ax.plot(midpoints, mean, color=color,
#                linestyle=linestyle, label='mean', **kwargs)
#        ax.set_xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
#        ax.set_ylabel(r"HOLE radius $R$ ($\AA$)")
#        if legend:
#            ax.legend(loc=legend_loc)
#        return ax

