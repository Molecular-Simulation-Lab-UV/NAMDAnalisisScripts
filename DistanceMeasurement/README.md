[![Prody](https://img.shields.io/badge/powered%20by-ProDy-green)](http://prody.csb.pitt.edu/index.html)
# Distance Measurement Scripts for NAMD trajectories

## distanceCalculation.py
WARNING: This script takes advantage of multicore processing due to `prody` and `numpy` configuration.
If used in don-elias, set the Open MP threads variable first to limit core usage: `export OMP_NUM_THREADS=6` should be enough to get good performance.

&nbsp;&nbsp; |- Running: `python distanceCalculation.py -i input_file [-d] dimension` <br />
&nbsp;&nbsp; |- Reminder/Help: `python distanceCalculation.py -h` <br />
&nbsp;&nbsp; |- Requires: Prody (can be installed via conda). I believe Prody requires MDAnalysis as well.

An example input file can be found together with the script. VMD style selections allow for spaces. File paths/names don't.
All paths, both via command line and in the input file, can be specified as relative or absolute (from the current folder, or from the root /)
