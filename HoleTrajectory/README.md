[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

# Batch HOLE analysis

## center.tcl
Small script to align trajectories. The user has to change the selection inside the script first. Then, usage is as follows (in don-elias), with vmd in the path: <br />
`vmd-1.9.3 -e center.tcl -args /path/to/psfFile.psf /path/to/refPDB.pdb /path/to/inputDCD.dcd/ /path/to/outputDCD.dcd`<br />
The input of `center.tcl` (inputDCD.dcd, third argument) should be a pre-processed .dcd file with only the pore/chain atoms to be analyzed, obtained from (in don-elias)

```
module load catdcd/5.1
catdcd -o outname.dcd -otype dcd -stride strideAmount -i atomIndexFile -first firstFrame -last lastFrame path/to/dcdFile.dcd
```
The `atomIndexFile` must be prepared beforehand, and should contain the pore/chain atom indexes only.
This output `outname.dcd` should be entered in the input file (inputTemplate.in) for the next script, under the "dcd" keyword.

## trajHOLE.py
This scripts leverages the MDAnalysis package for batch HOLE analysis of a .dcd, without having to convert every single frame to a .pdb file.
It's recommended that the .dcd trajectory be aligned to a reference (see `center.tcl` script described above),
and the pore to be studied should preferably be singled out, eliminating everything else appart from it from the trajectory.
Execution of the python script can be achieved by:
`python trajHOLE.py -i inputTemplate.in [-p /path/to/HOLEbinary]`
The input file for `trajHole.py` should contain:
&nbsp;&nbsp; |- `dcd`: Path to the pre-processed .dcd file. This file should preferably only contain the pore/chain along the trajectory. <br />
&nbsp;&nbsp; |- `psf`: Path to the corresponding .psf file. By that, we mean it should match the atoms contained in the .dcd file <br />
&nbsp;&nbsp; |- `bin`: Optional: A three-column argument. The first is the number of bins. The second and third are the lower and upper Z position for binning the pore. Defaults to [100, -20, 20] <br />
&nbsp;&nbsp; |- `sel`: Optional: A selection in (pseudo) VMD-like format. It's most useful when the .dcd trajectory is not correctly pre-processed to isolate the chain/pore. <br />
&nbsp;&nbsp; |- `out`: Optional: A path and/or name for the output file to be saved. Defaults to `outFile.dat` <br />
&nbsp;&nbsp; |- `keep`: Optional: A keyword to mark whether to keep or delete HOLE intermediate files. If present, files are kept. Remove the keyword to delete HOLE files when calculations are finished.
