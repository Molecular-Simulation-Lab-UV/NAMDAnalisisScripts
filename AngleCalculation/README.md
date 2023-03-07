[![Status](https://img.shields.io/badge/status-ready-green)]()
[![Prody](https://img.shields.io/badge/powered%20by-ProDy-9cf)](http://prody.csb.pitt.edu/index.html)
# Script for angle calculation.

## measureAngles.py
This script calculates the phi angle, the psi angle, or both, for a vmd-style selection specified in the input file.
BEFORE USING IT IN DON-ELIAS: `export OMP_NUM_THREADS = X`, where `X` is the number of cores to be used. Normally 6 would be enough to get results quickly, but it depends on the availability.

&nbsp;&nbsp; |- Running the script `python measureAngles.py -i input_file [-a angle]`. `angle` can be `phi`, `psi` or `both`.<br />
&nbsp;&nbsp; |- Reminder/Help: `python measureAngles.py -h`.<br />

If `out` is specified, this script will add the type of angle calculated as part of the output file's name. E.g.: out.out --> out.phi.out.
### Requirements
Prody. Prody also requires MD-Analysis as a core component.
