[![Status](https://img.shields.io/badge/status-testing-yellow)]()
[![Prody](https://img.shields.io/badge/powered%20by-ProDy-9cf)](http://prody.csb.pitt.edu/index.html)

# Script for angle calculation.

The script calculates the phi angle of the first selection residue only, and the psi angle of all residues in the selection.
Since the psi angle of residue j coincides with the phi angle of residue j+1, it's redundant to calculate both for intermediate residues.
