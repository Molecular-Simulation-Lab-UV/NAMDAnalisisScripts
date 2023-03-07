[![Powered by MDAnalysis](https://img.shields.io/badge/status-testing-yellow)](https://www.mdanalysis.org)

# Batch HOLE analysis

## center.tcl
Small script to align trajectories. The user has to change the selection inside the script first. Then, usage is as follows (in don-elias), with vmd in the path: <br />
```
vmd-1.9.3 -e center.tcl -args /path/to/psfFile.psf /path/to/refPDB.pdb /path/to/inputDCD.dcd/ /path/to/outputDCD.dcd
```
The input of `center.tcl` (inputDCD.dcd, third argument) should preferably be a pre-processed .dcd file with only the pore/chain atoms to be analyzed, obtained from (in don-elias, with catdcd loaded):

```
catdcd -o outname.dcd -otype dcd -stride strideAmount -i atomIndexFile -first firstFrame -last lastFrame path/to/dcdFile.dcd
```
The `atomIndexFile` must be prepared beforehand, and should contain the pore/chain atom indexes only.
This output `outname.dcd` should be entered in the input file (inputTemplate.in) for the next script, under the "dcd" keyword.

## trajHOLE.py
This script leverages the MDAnalysis package for batch HOLE analysis of a .dcd, without having to convert every single frame to a .pdb file.
It's recommended that the .dcd trajectory be aligned to a reference (see `center.tcl` script described above),
and the pore to be studied should preferably be singled out, eliminating everything else appart from it from the trajectory.
Execution of the python script can be achieved by running `python trajHOLE.py -i inputFile.in [-p /path/to/HOLEbinary]`<br/>
The input file for `trajHole.py` should contain:<br />
&nbsp;&nbsp; |- `dcd`: Path to the pre-processed .dcd file. This file should preferably only contain the pore/chain along the trajectory. <br />
&nbsp;&nbsp; |- `psf`: Path to the corresponding .psf file. By that, we mean it should match the atoms contained in the .dcd file <br />
&nbsp;&nbsp; |- `bin`: Optional: A three-column argument. The first is the number of bins. The second and third are the lower and upper Z position for binning the pore. Defaults to [100, -20, 20] <br />
&nbsp;&nbsp; |- `sel`: Optional: A selection in VMD-like format. It's most useful when the .dcd trajectory is not correctly pre-processed to isolate the chain/pore. Defaults to `protein`.<br />
&nbsp;&nbsp; |- `out`: Optional: A path and/or name for the output file to be saved. Defaults to `outFile.dat` <br />
&nbsp;&nbsp; |- `keep`: Optional: A keyword to mark whether to keep or delete HOLE intermediate files. If present, files are kept. Remove the keyword to delete HOLE files when calculations are finished. Defaults to `False` if not present.
