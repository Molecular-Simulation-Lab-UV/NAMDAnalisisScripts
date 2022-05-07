[![Prody](https://img.shields.io/badge/powered%20by-ProDy-9cf)](http://prody.csb.pitt.edu/index.html)
# Extract a selection from multiple .pdb files

## extractAsPDB.py
Allows reading of multiple .pdb files and creates a set of new .pdb files that contain only the selection entered as argument. Hence the name `extractAsPDB.py`
The new files are saved with, at most, the two first words of the selection appended to the original name, in the original directory.

Usage: `python extractAsPDB.py -p path/to/pdbDirectory [-s "VMD-style selection"]`

The only supported files format is `.pdb`, as it ignores all other extensions in the given directory.
