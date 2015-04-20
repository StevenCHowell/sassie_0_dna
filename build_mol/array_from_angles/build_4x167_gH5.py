#!/usr/bin/env python
#coding:utf-8
"""
  Author:  --<Steven C. Howell>
  Purpose: align monomer from 12mer with linker histone to c11
  Created: 08/06/2014
"""

import sassie.sasmol.sasmol as sasmol
import numpy as np
import time
import string, os, locale, sys, random
import x_dna.util.tetramer_angles as ta
import sassie_1_na.util.geometry as geometry

def renumber_index_resid(index, resid):
    '''
    generate a corrected index and resid list
    '''
    number=[]
    resid_array=[] ; count=1
    for i in xrange(len(resid)):
        this_resid = resid[i]
        if(i==0):
            last_resid = this_resid
        else:
            if(this_resid != last_resid):	
                count += 1
        resid_array.append(count)	
        last_resid = this_resid
        number.append(i+1)

    return index, resid


def combine_pdbs(all_pdbs, out_pdb=None):
    '''
    given a list of pdb files, this will combine them into one pdb
    
    inputs:
        all_pdbs - list of pdb file names
        out_pdb - optional file name to save the combined pdbs to
        
    outputs:
        combined_mol - the combined sasmol object
        
    see also:
        combine_sasmols
    '''
    all_mols = []
    for (i, pdb) in enumerate(all_pdbs):
        mol = sasmol.SasMol(0)
        mol.read_pdb(pdb)
        all_mols.append(mol)
        
    combined_mol = combine_sasmols(all_mols)
    if out_pdb: combined_mol.write_pdb(out_pdb, 0, 'w')
    return combined_mol
    
def combine_sasmols(all_mols, combine_segnames=None):
    '''
    given a list of sasmol objects, this will combine them into one sasmol object
    
    inputs:
        all_mols - list of sasmol objects
        
    outputs:
        combined_mol - the combined sasmol object
        
    see also:
        combine_pdbs
    '''
    combined_mol = sasmol.SasMol(0)
    for (i, mol) in enumerate(all_mols):
        i_str = str(i + 1)
        if i == 0:
            combined_mol.setAtom(mol.atom())
            combined_mol.setIndex(mol.index())
            combined_mol.setName(mol.name())
            combined_mol.setLoc(mol.loc())
            combined_mol.setResname(mol.resname())
            combined_mol.setChain(mol.chain())
            combined_mol.setResid(mol.resid())
            combined_mol.setRescode(mol.rescode())
            combined_mol.setCoor(mol.coor())
            combined_mol.setOccupancy(mol.occupancy())
            combined_mol.setBeta(mol.beta())
            combined_mol.setSegname([i_str + name for name in mol.segname()])
            combined_mol.setElement(mol.element())
            combined_mol.setCharge(mol.charge())
        else:
            index = np.concatenate((combined_mol.index(), combined_mol.index()[-1] + mol.index()))
            combined_mol.setIndex(index)
            combined_mol.setAtom(combined_mol.atom() + mol.atom())
            combined_mol.setName(combined_mol.name() + mol.name())
            combined_mol.setLoc(combined_mol.loc() + mol.loc())
            combined_mol.setResname(combined_mol.resname() + mol.resname())
            combined_mol.setChain(combined_mol.chain() + mol.chain())
            resid = np.concatenate((combined_mol.resid(), combined_mol.resid()[-1] + mol.resid()))
            combined_mol.setResid(resid)
            combined_mol.setRescode(combined_mol.rescode() + mol.rescode())
            coor = np.concatenate((combined_mol.coor(), mol.coor()), axis=1)
            combined_mol.setCoor(coor)
            combined_mol.setOccupancy(combined_mol.occupancy() + mol.occupancy())
            combined_mol.setBeta(combined_mol.beta() + mol.beta())
            combined_mol.setSegname(combined_mol.segname() + [i_str + name for name in mol.segname()])
            combined_mol.setElement(combined_mol.element() + mol.element())
            combined_mol.setCharge(combined_mol.charge() + mol.charge())

    segname = [name.replace('DA', 'ADE') for name in combined_mol.resname()]
    
    if combine_segnames:
        for segname in combine_segnames:
            tmp_segname = combined_mol.segname()
            for j in xrange(i+1):
                j_str = str(j+1)
                tmp_segname = [name.replace(j_str + segname, segname) for name in tmp_segname]
            
            combined_mol.setSegnames(tmp_segnames)
            error, mask = combined_mol.get_subset_mask("segname[i] = '%s'" % segname)
            index, resid = renumber_index_resid(index, resid)
            combined_mol.setIndex(index)
            combined_mol.setResid(resid)    

    return combined_mol

def align_gH5_to_c11():
    # load in the C11 pdb   
    c11 = sasmol.SasMol(0)
    c11.read_pdb('c11.pdb')

    # create the filter to select the C11 DNA to align to
    c11_filter = []
    c11_filter.append("((name[i] == 'P') and (segname[i] == 'DNA1') and (resid[i] <= 166) and (resid[i] >=  26))")  # NCP-1
    c11_filter.append("((name[i] == 'P') and (segname[i] == 'DNA1') and (resid[i] <= 333) and (resid[i] >= 193))")  # NCP-2
    c11_filter.append("((name[i] == 'P') and (segname[i] == 'DNA1') and (resid[i] <= 499) and (resid[i] >= 359))")  # NCP-3
    c11_filter.append("((name[i] == 'P') and (segname[i] == 'DNA1') and (resid[i] <= 667) and (resid[i] >= 527))")  # NCP-4

    all_ncps = []

    for i in xrange(len(c11_filter)):
        # load in the nucleosome with the linker protein
        linker = sasmol.SasMol(0)
        linker.read_pdb('gH5_NCP.pdb')

        linker.center(0)
        # create the filter to select the coordinates for aligning
        linker_filter = "((moltype[i] == 'dna' or moltype[i] == 'rna') and name[i] == 'P' and chain[i] == 'I' and resid[i] >= 14 and resid[i] <= 154)"
        error, linker_mask = linker.get_subset_mask(linker_filter)
        s_link = sum(linker_mask)
        print 'sum(linker_mask) =', s_link 

        # load the alignment coordinates into its own sosmol object
        linker_sub = sasmol.SasMol(0)
        error = linker.copy_molecule_using_mask(linker_sub, linker_mask, 0)
        linker.center(0)  #m2

        # (sub_m2) get the com and coordinates for each protein to be used for alignment
        com_sub_2 = linker_sub.calccom(0)
        linker_sub.center(0)
        coor_sub_2 = linker_sub.coor()[0]

        linker.center(0)  # m2: linker
        coor_linker = linker.coor()

        print "aligning to NCP", i

        error, c11_mask = c11.get_subset_mask(c11_filter[i])
        s_c11 = sum(c11_mask)
        # print 'sum(c11_mask) =',  s_c11

        # copy the protein to its own sasmol object
        c11_sub = sasmol.SasMol(0)
        error = c11.copy_molecule_using_mask(c11_sub, c11_mask, 0)

        com_sub_1 = c11_sub.calccom(0)
        c11_sub.center(0)
        coor_sub_1 = c11_sub.coor()[0]

        assert s_link == s_c11, "c11 has %d atoms but H5 has %d atoms" % (s_c11, s_link)
        linker.align(0, coor_sub_2, com_sub_2, coor_sub_1, com_sub_1)
        # sub_1: the molecule to align too
        # sub_2: the molecule to bp aligned

        # write the output to a new file
        ncp_filename = "NCP%d.pdb" % (i + 1)
        linker.write_pdb(ncp_filename, 0, 'w')

        all_ncps.append(linker)

    # ncp_array = combine_sasmols(all_ncps, combine_segnames=['I', 'J'])
    ncp_array = combine_sasmols(all_ncps)
    resname = [name.replace('DA', 'ADE').replace('DC', 'CYT').replace('DG', 'GUA').replace('DT', 'THY') for 
               name in ncp_array.resname()]
    ncp_array.setResname(resname)
    pdb = '4x167_gH5.pdb'
    ncp_array.write_pdb(pdb, 0, 'w') 
    
    bps = np.array([np.linspace(0, 167, 168), np.linspace(168, 1, 168)]).T
    ncp_dna_resids = [bps[[14, 154]]]*4
    dna_ids = [['1I', '1J'], ['2I', '2J'], ['3I', '3J'], ['4I', '4J']]
    ncp_dyad_resids = [bps[84], bps[84], bps[84], bps[84]]
    all_ncp_plot_vars, all_ncp_axes, all_ncp_origins = ta.get_tetramer_axes(pdb, ncp_dna_resids, dna_ids, ncp_dyad_resids)
    geometry.show_ncps(all_ncp_plot_vars)

    print '\m/ >.< \m/'
    '''
    tic = time.time()
    dna_filter = "((moltype[i] == 'dna') or (moltype[i] == 'rna'))"
    error, dna_mask = c11.get_subset_mask(dna_filter)
    dna = sasmol.SasMol(0)
    error = c11.copy_molecule_using_mask(dna, dna_mask, 0)
    toc_dna = time.time() - tic

    tic = time.time()
    pro_filter = "(name[i] == 'CA')"
    error, pro_mask = c11.get_subset_mask(pro_filter)
    pro = sasmol.SasMol(0)
    error = c11.copy_molecule_using_mask(pro, pro_mask, 0)
    toc_pro = time.time() - tic

    tic = time.time()
    not_dna = sasmol.SasMol(0)
    not_dna_mask = ((dna_mask-1) * -1 )
    error = c11.copy_molecule_using_mask(not_dna, not_dna_mask, 0)   
    pro_filter2 = "(name[i] == 'CA')"
    error, pro_mask2 = not_dna.get_subset_mask(pro_filter2)
    pro2 = sasmol.SasMol(0)
    error = not_dna.copy_molecule_using_mask(pro2, pro_mask2, 0)
    toc_not_dna = time.time() - tic



    print 'dna:', toc_dna
    print 'pro:', toc_pro
    print 'not_dna:', toc_not_dna
    print 'sum dna_mask:', sum(dna_mask)
    print 'sum pro_mask:', sum(pro_mask)
    print 'sum pro_mask2:', sum(pro_mask2)
    # test = (pro_mask_i - dna_mask) / 2
    # print 'max(test):', max(test)
    '''
    
    
if __name__ == '__main__':
    n_ncp = 4
    align_gH5_to_c11()
    