'''
This file is meant to serve as a record of the generation of the dimer.pdb file
rather than a python library to be improted.  It also serves as an example of 
how to use several of the libraries called within.
'''
import sassie.sasmol.sasmol as sasmol
import x_dna.util.basis_to_python as basis_to_python
import numpy as np
import x_dna.drivers.myAlign as align
    
class inputs():
    def __init__(self, parent = None):
        pass

def main():
    ''' 
    this did not work well, there histones to not match well 
    in the NCP so there was overlap btwn the aligned versions
    '''
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
    dimer_file = '../150205dimer.pdb'

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
    in_vars.pdb = '../dna2_right.pdb'
    in_vars.segnames = ['dummy','DNA2']
    pdb2psfgen.main(in_vars)

    in_vars.pdb = '../dna1_right.pdb'
    in_vars.segnames = ['DNA1','dummy']
    pdb2psfgen.main(in_vars)
    
def replace_N_atoms():
    '''
    swap N9-N1 atoms in rename DNA residues
    GUA or ADE N1 -> N9
    CYT or THY N1 -> N9    
    '''
    dna1_file = '../dna1_right_seq.pdb'
    dna1_out = '../dna1_right.pdb'
    dna1 = sasmol.SasMol(0)
    dna1.read_pdb(dna1_file)


    dna2_file = '../dna2_right_seq.pdb'
    dna2_out = '../dna2_right.pdb'
    dna2 = sasmol.SasMol(0)
    dna2.read_pdb(dna2_file)    
    
    '''
    Count the number of resids for each resid
    Loop over each resid with <13 atoms 
    .  Use enumerate to get the indices for the atoms the residue 
    .  Iterate over those indices to find the N1 or N9 atom
    .  Depending on which base type it is, replace the atom name 
    Store the names in the Sasmol object 
    Save the pdb

    Loop over every atom,
    add 1 to the number of some in that residue (can use count instead) 
    Store the index for the n9 and n1 atoms 
    Then loop over just the group s to for n atoms 
    '''

    replace_n1_n9(dna1, dna1_out)
    replace_n1_n9(dna2, dna2_out)

    return

def replace_n1_n9(dna_mol, dna_out):
    import pandas as pd

    pyrimidines = ['CYT', 'THY']
    purines = ['ADE', 'GUA']
    n_atoms = [(dna_mol.resid()==resid).sum() for resid in dna_mol.resids()]
    i_res_min = []
    i_res_max = []
    for j in dna_mol.resids():
        i_res = [i for i, resid in enumerate(dna_mol.resid()) if resid == j]
        i_res_min.append(min(i_res))
        i_res_max.append(max(i_res))
    
    resname = [dna_mol.resname()[i] for i in i_res_min]
    zeros = np.zeros(len(n_atoms),dtype=int)
    dna_dict = {'n_atoms': n_atoms, 'resid': dna_mol.resids(), 'resname': resname,
                 'i_N9': zeros, 'i_N1': zeros, 'i_min': i_res_min, 
                 'i_max': i_res_max, 'i_replace': zeros}
    headers = ['resid', 'n_atoms', 'resname', 'i_min',
               'i_max', 'i_N9', 'i_N1', 'i_replace']
    dna_frame = pd.DataFrame(dna_dict, columns=headers)
    
    N_frame = dna_frame[dna_frame['n_atoms'] < 13].reset_index()
    names = dna_mol.name()
    before = dna_mol.name()
    for i in xrange(len(N_frame)):
        i_min = N_frame['i_min'][i]
        i_max = N_frame['i_max'][i] + 1
        try:
            N_frame['i_N1'][i] = names[i_min:i_max].index('N1')
            N_frame['i_N9'][i] = -1
            N_frame['i_replace'][i] = N_frame['i_N1'][i] + N_frame['i_min'][i]
        except:
            N_frame['i_N9'][i] = names[i_min:i_max].index('N9')
            N_frame['i_N1'][i] = -1
            N_frame['i_replace'][i] = N_frame['i_N9'][i] + N_frame['i_min'][i]

            
        if N_frame['resname'][i] in pyrimidines and N_frame['i_N9'][i] >= 0:
            names[N_frame['i_replace'][i]] = 'N1'
            # print 'did switch %s, N1(%d), N9(%d)' % (N_frame['resname'][i],
                                                         # N_frame['i_N1'][i],
                                                         # N_frame['i_N9'][i])

        elif N_frame['resname'][i] in purines and N_frame['i_N1'][i] >= 0:
            names[N_frame['i_replace'][i]] = 'N9'
            # print 'did switch %s, N1(%d), N9(%d)' % (N_frame['resname'][i],
                                                         # N_frame['i_N1'][i],
                                                         # N_frame['i_N9'][i])
        
        else:
            print 'did not switch %s, N1(%d), N9(%d)' % (N_frame['resname'][i],
                                                         N_frame['i_N1'][i],
                                                         N_frame['i_N9'][i])

    dna_mol.setName(names)
    after = dna_mol.name()
    if dna_mol.write_pdb(dna_out, 0, 'w'):
        print 'sucessfully wrote %s' % dna_out
    
if __name__ == "__main__":
    
    # align_1kx5()
    # separate_1kx5()
    # verify_dna_sequence()
    # replace_dna_sequence()
    generate_psfgen_patches()
    # main()
    # replace_N_atoms()
    
    print 'done'
