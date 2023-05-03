[![Status](https://img.shields.io/badge/status-working-green)]()
[![Prody](https://img.shields.io/badge/powered%20by-ProDy-9cf)](http://prody.csb.pitt.edu/index.html)

# Permeation Events count

## Permeation events
To count permeation events, we consider a two cylinders as the definition of the entry points, both upper and lower, through which the atom or molecule of interest permeates. The system is split into 5 zones, namely a lower-outer volume (flag 1), the upper outer volume (flag 5), the two cylinders previously mentioned (flags 2 and 4 for the lower and upper one, respectively), and a middle slab that is bounded only in z by the inner limits of both cylinders (flag 3). So, we define permeation events as the transit of an atom or molecule from one outer volume to the other, passing through both cylinders and the central slab. The cylinders' long axes are always parallel to the simulations' Z axis, and the radius is on a plane parallel to that of the membrane (XY plane). <br />

The script makes extensive use of ProDy, harnessing its coordinate obtainment methods, as well as its wrapping and fitting methods. The algorithm is rather simple, tagging water molecules each timestep, updating a flags "history" (array `flagsOld`) each time a new flag is found for an element of the selection, and resetting the array if the element is back in bulk. If the element being studied travels in such a way that all flags are present in a specific order in the array, a permeation event is counted.

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
&nbsp;&nbsp; |- `delta`: \[Optional] Height (z-length) of each cylinder defining the entry/exit points of the pore. Default is 3\[Å] <br />
&nbsp;&nbsp; |- `rad`: Radius of the bounding cylinder, measured in \[Å] <br />
&nbsp;&nbsp; |- `ref`: VMD-like selection of the elements used for alignment between each frame in the trajectory and the reference structure. <br/>
&nbsp;&nbsp; |- `sel`: VMD-like selection of the atoms or molecules whose permeation will be counted. E.g: `resname CLA` .<br/>
&nbsp;&nbsp; |- `out`: \[Optional] Path and name of the output file. Default is `outFile.out`.
