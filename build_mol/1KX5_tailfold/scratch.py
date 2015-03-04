import numpy as np
import sassie.sasmol.sasmol as sasmol
import sassie_1_na.util.basis_to_python as basis_to_python
import x_dna.build_mol as build_mol

ncp = sasmol.SasMol(0)
ncp.read_pdb('1KX5tailfold_167bp.pdb')



# protein chains
# H3: A, E # E will conflict with the gH5 protein, try aligning A instead
# H4: B, F
# H2A: C, G # Missing an A residue at the end
# H2B: D, H
#
# dna chains
# I, J, K, L, M, N

# todo:
# 1. align the 1KX5 proteins to the correct positions
# X. re-confirm the sequence files from before, rename their chain names to:
#    A: H2A
#    B: H2B
#    T: H3
#    F: H4
# 3. create pdb and psf files for starting structures, name the segments:
#    1H2A - 8H2A
#    1H2B - 8H2B
#    1H3  - 8H3
#    1H4  - 8H4
#      - dimer 2x167
#      - trimer 3x167
#      - tetramer 4x167
#      - gH5 tetramer 4x167
# 4. create driver script for each then run making the following variations:
#    a. amount of unwrapping
#    b. twist angle (add in the twist energy)
#    c. bend angle
