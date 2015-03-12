'''
This file is meant to serve as a record of the generation of the dimer.pdb file
rather than a python library to be improted.  It also serves as an example of 
how to use several of the libraries called within.
'''
import sassie.sasmol.sasmol as sasmol
import x_dna.util.basis_to_python as basis_to_python

class inputs():
    def __init__(self, parent = None):
        pass

def main():
    import x_dna.drivers.myAlign as align
    in_vars = inputs()
    
    goals = ['histone_A0.pdb', 'histone_B0.pdb', 'histone_C0.pdb', 'histone_D0.pdb', 
             'histone_E0.pdb', 'histone_F0.pdb', 'histone_G0.pdb', 'histone_H0.pdb', 
             'histone_A1.pdb', 'histone_B1.pdb', 'histone_C1.pdb', 'histone_D1.pdb', 
             'histone_E1.pdb', 'histone_F1.pdb', 'histone_G1.pdb', 'histone_H1.pdb']
    move = ['1KX5tailfold_167bp_chain_A1.pdb', '1KX5tailfold_167bp_chain_B1.pdb', 
            '1KX5tailfold_167bp_chain_C1_trun.pdb', '1KX5tailfold_167bp_chain_D1.pdb', 
            '1KX5tailfold_167bp_chain_A1.pdb', '1KX5tailfold_167bp_chain_F1.pdb', 
            '1KX5tailfold_167bp_chain_G1_trun.pdb', '1KX5tailfold_167bp_chain_H1.pdb']
    out = ['1H3.pdb', '1H4.pdb', '1H2A.pdb', '1H2B.pdb',
           '2H3.pdb', '2H4.pdb', '2H2A.pdb', '2H2B.pdb',
           '3H3.pdb', '3H4.pdb', '3H2A.pdb', '3H2B.pdb',
           '4H3.pdb', '4H4.pdb', '4H2A.pdb', '4H2B.pdb']
    move_seg = ['A', 'B', 'C', 'D', 'A', 'F', 'G', 'H']
    goal_seg = ['A0', 'B0', 'C0', 'D0', 
                'E0', 'F0', 'G0', 'H0', 
                'A1', 'B1', 'C1', 'D1', 
                'E1', 'F1', 'G1', 'H1']
    in_vars.move_seg_or_ch = 'segname'
    in_vars.goal_seg_or_ch = 'segname'
    in_vars.path = './'
    match_res_min = [50, 20, 01, 50]
    match_res_max = [70, 40, 50, 124]
    
    for i in xrange(len(goals)):
        in_vars.goal = goals[i]
        in_vars.move = in_vars.ref = move[i%8]
        in_vars.out = out[i]
        in_vars.move_seg_chain = move_seg[i%8]
        in_vars.goal_seg_chain = goal_seg[i]
        in_vars.min = match_res_min[i%4]
        in_vars.max = match_res_max[i%4]
    
    align.align(in_vars)

def mv_A_to_E():
    '''
    this worked well to move A to E (no apparent overlap)
    '''
    in_vars = inputs()
    in_vars.goal = '1KX5tailfold_167bp_chain_E1.pdb'
    in_vars.ref  = '1KX5tailfold_167bp_chain_A1.pdb'
    in_vars.move = '1KX5tailfold_167bp_chain_A1.pdb'    
    in_vars.out  = '1KX5tailfold_A2E_2.pdb'
    in_vars.move_seg_chain = 'A'
    in_vars.goal_seg_chain = 'E'
    in_vars.path = './'
    in_vars.move_seg_or_ch = 'segname'
    in_vars.goal_seg_or_ch = 'segname'    
    in_vars.min = 38
    in_vars.max = 50
    
    align.align(in_vars)
    

def replace_tail(sasmol, basis_to_python):
    # mv_A_to_E()
    
    mono_file = '../../1KX5_tailfold/1KX5tailfold_167bp.pdb'
    ncp = sasmol.SasMol(0)
    ncp.read_pdb(mono_file)
    replace_basis = basis_to_python.parse_basis('chain E and resid < 40')
    error, replace_mask = ncp.get_subset_mask(replace_basis)
    print sum(replace_mask)
    
    histone_file = '1KX5tailfold_A2E_2.pdb'
    new_histone = sasmol.SasMol(0)
    new_histone.read_pdb(histone_file)
    print new_histone.coor().shape
    part_histone = sasmol.SasMol(0)
    
    part_basis = basis_to_python.parse_basis('resid < 40')
    error, part_mask = new_histone.get_subset_mask(part_basis)
    new_histone.copy_molecule_using_mask(part_histone, part_mask, 0)
    
    ncp.set_coor_using_mask(part_histone, 0, replace_mask)
    ncp.write_pdb('1KX5tailfold_fxd2.pdb', 0, 'w')
    
def align_1kx5():
    '''
    align the 1kx5 to the dimer then save pdbs for psfgen
    '''
    mono_file = '1KX5tailfold_fxd.pdb'
    dimer_file = '../dimer.pdb'

    in_vars = inputs()
    in_vars.goal = dimer_file
    in_vars.ref = in_vars.move = mono_file
    in_vars.path = './'
    in_vars.move_filter = '((chain[i] == "I") and (name[i] == "C1\'"))'

    # NCP 1
    in_vars.out = '1KX5tailfold_dimer_ncp1.pdb'
    in_vars.goal_filter = '((segname[i] == "DNA1") and (name[i] == "C1\'") and (resid[i] > 14) and (resid[i] < 162))'
    align.align(in_vars)

    # NCP 2    
    in_vars.out = '1KX5tailfold_dimer_ncp2.pdb'
    in_vars.goal_filter = '((segname[i] == "DNA1") and (name[i] == "C1\'") and (resid[i] > 181) and (resid[i] < 329))'
    align.align(in_vars)
    
def separate_1kx5():
    '''
    create separate pdb files of all the histones
    '''
    import x_dna.build_mol.chain_get_pdb_and_seq as get_pdb
    in_vars = inputs()
    in_vars.chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    in_vars.chains = ['C', 'G']

    in_vars.pdb = 'dimer_ncp1/1KX5tailfold_dimer_ncp1.pdb'
    get_pdb.main(in_vars)

    in_vars.pdb = 'dimer_ncp2/1KX5tailfold_dimer_ncp2.pdb'
    get_pdb.main(in_vars)
    
def verify_dna_sequence():
    '''
    checking the dna sequence of the dna input files
    '''
    import x_dna.build_mol.seg_get_pdb_and_seq as get_seq

    in_vars = inputs()

    in_vars.segnames = ['DNA1']
    in_vars.pdb = '../dna1_right_seq.pdb'
    get_seq.main(in_vars)

    in_vars.segnames = ['DNA2']
    in_vars.pdb = '../dna2_right_seq.pdb'
    get_seq.main(in_vars)
    # this one was completely wrong !!!
    
def replace_dna_sequence():
    '''
    replace the incorrect dna sequence is the DNA2 segment
    '''
    import x_dna.build_mol.replace_sequence as replace_seq
    replace_seq.main()
    
def generate_psfgen_patches():
    '''
    generate the psfgen patches for the dimer dna chains
    '''
    import x_dna.build_mol.seg_pdb2psfgen as pdb2psfgen
    in_vars = inputs()
    in_vars.pdb = '../dna2_bb_right-seq.pdb'
    in_vars.segnames = ['dummy','DNA2']
    pdb2psfgen.main(in_vars)

    in_vars.pdb = '../dna1_bb_right-seq.pdb'
    in_vars.segnames = ['DNA1','dummy']
    pdb2psfgen.main(in_vars)
    
def replace_N_atoms():
    '''
    swap N9-N1 atoms in rename DNA residues
    GUA or ADE N1 -> N9
    CYT or THY N1 -> N9    
    '''
    dna1_file = '../dna1_right_seq.pdb'
    dna2_file = '../dna2_right_seq.pdb'
    dna1 = sasmol.SasMol(0)
    dna2 = sasmol.SasMol(0)
    dna1.read_pdb(dna1_file)
    dna2.read_pdb(dna2_file)    
    
    natoms = []
    pyrimidines = ['CYT', 'THY']
    purines = ['ADE', 'GUA']
    n1_basis = basis_to_python.parse_basis('name N1')
    residue = sasmol.SasMol(0)
    for resid in dna1.resids():
        res_basis = basis_to_python.parse_basis('resid %d' % resid)
        error, res_mask = dna1.get_subset_mask(res_basis)
        natoms.append(res_mask.sum())
        if res_mask.sum() < 20:
            dna1.copy_molecule_using_mask(residue, res_mask, 0)
            if residue.resnames()[0] in purines:
                names = residue.name()
                for (i, name) in enumerate(names):
                    if name == 'N1':
                        
                        names[i] = 'N9'
                residue.setName(names)
            elif residue.resnames()[0] in pyrimidines:
                names = residue.name()
                for (i, name) in enumerate(names):
                    if name == 'N9':
                        names[i] = 'N1'
                residue.setName(names)
            dna1.set_

    N_atoms = []
    current_resid = 0
    natoms_in_current_res = 0
    np.zeros(len(r))  #populate this <---
    natoms_in_res = []
    for i in dna1.index()[:40]:
        if dna1.resid()[i] == current_resid:
            natoms_in_current_res += 1
        else:
            natoms_in_res.append(n_atoms_in_current_res)
            current_resid = dna1.resid()[i]
            natoms_in_current_res = 1
            
        if dna1.name()[i] in ['N1', 'N9']:
            N_atoms.append(i)
        
            

                
    dna1_out = '../dna1_right.pdb'    
    dna2_out = '../dna2_right.pdb'    
    dna1.write_pdb(dna1_out, 0, 'w')
    dna2.write_pdb(dna2_out, 0, 'w')
    return
    
if __name__ == "__main__":
    
    # align_1kx5()
    # separate_1kx5()
    # verify_dna_sequence()
    # replace_dna_sequence()
    # generate_psfgen_patches()
    # main()
    replace_N_atoms()
    
    print 'done'