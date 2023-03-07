![](https://img.shields.io/badge/status-testing-yellow)
[![Prody](https://img.shields.io/badge/powered%20by-ProDy-9cf)](http://prody.csb.pitt.edu/index.html)

# Permeation Events count

## Permeation events
To count permeation events, we consider a cylinder as the definition of the pore through which the atom or molecule of interest permeates. So, we define permeation events as the transit of an atom or molecule from one side of that cylinder to the other, passing through said cylinder. Its long axis is always parallel to the simulations' Z axis, and the radius is on a plane parallel to that of the membrane (XY plane). <br />

## Running the program
REMEMBER TO LIMIT YOUR VISIBLE CPUs!<br />
`export OMP_NUM_THREADS=6` should be enough.

The script wraps trajectories, centers them, and analyzes them.<br/>
Execution of the python script can be achieved by running `python permeationEvents.py -i inputFile.in`<br/>
The input file for `permeationEvents.py` should contain:<br />
&nbsp;&nbsp; |- `dcd`: Path to the `.dcd` file. WARNING: This has not been tested with an already wrapped trajectory. <br />
&nbsp;&nbsp; |- `pdb`: Path to the corresponding .pdb file. <br />
&nbsp;&nbsp; |- `upperZ`: Upper Z coordinate of the pore's bounding cylinder in the reference pdb structure. <br />
&nbsp;&nbsp; |- `lowerZ`: Lower Z coordinate of the pore's bounding cylinder in the reference pdb structure. <br />
&nbsp;&nbsp; |- `rad`: Radius of the bounding cylinder, measured in \[Å] <br />
&nbsp;&nbsp; |- `ref`: VMD-like selection of the elements used for alignment between each frame in the trajectory and the reference structure. <br/>
&nbsp;&nbsp; |- `sel`: VMD-like selection of the atoms or molecules whose permeation will be counted. E.g: `resname CLA` .<br/>
&nbsp;&nbsp; |- `out`: \[Optional] Path and name of the output file. Default is `outFile.out`.
