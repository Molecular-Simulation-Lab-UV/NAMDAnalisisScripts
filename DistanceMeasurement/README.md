[![Prody](https://img.shields.io/badge/powered%20by-ProDy-9cf)](http://prody.csb.pitt.edu/index.html)
# Distance Measurement Scripts for NAMD trajectories

## distanceCalculation.py
WARNING: This script takes advantage of multicore processing due to `prody` and `numpy` configuration.
If used in don-elias, set the Open MP threads variable first to limit core usage: `export OMP_NUM_THREADS=6` should be enough to get good performance.

This scripts calculates the one of 3 options: The Z distance (-d 1), the radial distance (-d 2), or the 3D distance (default, -d 3). It outputs a file with the number of selections minus 1 (N - 1) columns of distances. The first column is the frame index. The selection labels are used as headers for each column, which might make the output file's first line a bit clumpy. Distances are calculated between the first selection `sel`, from top to bottom, and all others. So, this is a 1 to N calculation.

&nbsp;&nbsp; |- Running: `python distanceCalculation.py -i input_file [-d dimension]` <br />
&nbsp;&nbsp; |- Reminder/Help: `python distanceCalculation.py -h` <br />
&nbsp;&nbsp; |- Requires: Prody (can be installed via conda). I believe Prody requires MDAnalysis as well.

An example input file can be found together with the script. VMD style selections allow for spaces. File paths/names don't.
All paths, both via command line and in the input file, can be specified as relative or absolute (from the current folder, or from the root /)<br />
Input file arguments: <br />
&nbsp;&nbsp; |- `pdb`: Path to the reference .pdb file. This is used to align the trajectory frames.<br />
&nbsp;&nbsp; |- `ref`: VMD-style selection to be used for aligning the frames to the reference pdb. Used without "". For example: `chain D and backbone`<br />
&nbsp;&nbsp; |- `dcd(1/2/3/etc)`: Paths to .dcd files. The index in the keyword itself is not important, as the files will be sorted numerically within the script. Minimum is one .dcd file path. <br />
&nbsp;&nbsp; |- `sel(1/2/3/etc)`: VMD-like format selections. Again, the index (1/2/3/etc) in itself is not important. What matters is the order, top to bottom, of the keywords in the file. Minimum amount: 2 selections. To be written without "". For example: `resid 180 and name CA`.<br/>
&nbsp;&nbsp; |- `out`: Optional: Path and name for the output file. Defaults to `outFile.dat`
