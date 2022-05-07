# Distance Measurement Scripts for NAMD trajectories

## distanceCalculation.py
Usage: `python distanceCalculation.py -i input_file [-d] dimension`\n
Reminder/Help: `python distanceCalculation.py -h`
Requires: Prody (can be installed via conda). I believe Prody requires MDAnalysis as well.

An example input file can be found together with the script. VMD style selections allow for spaces. File paths/names don't.
All paths, both via command line and in the input file, can be specified as relative or absolute (from the current folder, or from the root /)
