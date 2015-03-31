'''
generate the psfgen patches for the dimer dna chains
'''
class inputs():
    def __init__(self, parent = None):
        pass

import x_dna.build_mol.seg_pdb2psfgen as pdb2psfgen
in_vars = inputs()

in_vars.pdb = 'dna_pdb/w601dna.pdb'
in_vars.segnames = ['DNA1','DNA2']
pdb2psfgen.main(in_vars)
