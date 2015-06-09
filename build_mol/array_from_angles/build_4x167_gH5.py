#!/usr/bin/env python
#coding:utf-8
"""
  Author:  Steven C. Howell
  Purpose: manually construct gH5 nucleosome arrays
  Created: 20 April 2015
  $Id$
00000000011111111112222222222333333333344444444445555555555666666666677777777778
12345678901234567890123456789012345678901234567890123456789012345678901234567890
"""

import sassie.sasmol.sasmol as sasmol
import numpy as np
# import time
# import string, os, locale, sys, random
import x_dna.util.tetramer_angles as ta
import x_dna.util.geometry as geometry
import x_dna.drivers.myAlign as align
import x_dna.util.basis_to_python as b2p

class inputs():
    def __init__(self, parent = None):
        pass

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

# m.setName([name.replace("*","'") for name in m.name()])
# def replace_prime(names):
    # last_name = ''
    # for (j, name) in enumerate(names):
        # if name == last_name: 
            # pass
        # else:
            # last_name = name
            # new_name = name.replace("*","'")
        # segnames[j] = new_name
        
def increment_ncp_segnames(segnames, duplicate_segnames, i):
    if len(duplicate_segnames):
        last_name = ''
        for (j, name) in enumerate(segnames):
            if name in duplicate_segnames:
                if name == last_name: 
                    pass
                else:
                    last_name = name
                    try: 
                        n = int(name[0])
                        new_name = str(n + 2*i) + name[1:]
                    except:
                        n = int(name[-1])
                        new_name = name[:-1] + str(n + 2*i)
                segnames[j] = new_name
    
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
            combined_mol.setSegname(mol.segname())
            combined_mol.setElement(mol.element())
            combined_mol.setCharge(mol.charge())
            combined_mol.setMoltype([moltype.replace('rna','dna') for moltype in mol.moltype()])
            combined_mol.setSegnames(list(set(combined_mol.segname())))
        else:
            index = np.concatenate((combined_mol.index(), combined_mol.index()[-1] + mol.index()))
            combined_mol.setIndex(index)
            combined_mol.setAtom(combined_mol.atom() + mol.atom())
            combined_mol.setName(combined_mol.name() + mol.name())
            combined_mol.setLoc(combined_mol.loc() + mol.loc())
            combined_mol.setResname(combined_mol.resname() + mol.resname())
            combined_mol.setChain(combined_mol.chain() + mol.chain())
            resid = np.concatenate((combined_mol.resid(), mol.resid()))
            combined_mol.setResid(resid)
            combined_mol.setRescode(combined_mol.rescode() + mol.rescode())
            coor = np.concatenate((combined_mol.coor(), mol.coor()), axis=1)
            combined_mol.setCoor(coor)
            combined_mol.setOccupancy(combined_mol.occupancy() + mol.occupancy())
            combined_mol.setBeta(combined_mol.beta() + mol.beta())
            duplicate_segnames = []
            for segname in mol.segnames():
                if segname in combined_mol.segnames():
                    duplicate_segnames.append(segname)
            increment_ncp_segnames(mol.segname(), duplicate_segnames, i)
            combined_mol.setSegname(combined_mol.segname() + mol.segname())
            combined_mol.setElement(combined_mol.element() + mol.element())
            combined_mol.setCharge(combined_mol.charge() + mol.charge())
            combined_mol.setMoltype(combined_mol.moltype() + 
                                    [moltype.replace('rna','dna') for moltype in mol.moltype()])
            combined_mol.setSegnames(list(set(combined_mol.segname())))
    
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
    pdb = '4x167_gH5.pdb'

    try:
        ncp_array = sasmol.SasMol(0)
        ncp_array.read_pdb(pdb)
    except:
        # load in the C11 pdb   
        c11 = sasmol.SasMol(0)
        c11.read_pdb('c11_r.pdb')
    
        # c11.read_pdb('c11.pdb')
        # # re-orient the array
        # coor = c11.coor()
        # coor[0] = geometry.transform_coor(coor[0], np.array([1, 0, 0]), np.array([0, 0, 0]))
        # c11.setCoor(coor)
        # c11.write_pdb('c11_r.pdb', 0, 'w')
        
        # create the filter to select the C11 DNA to align to
        c11_filter = []
        c11_filter.append("((name[i] == 'P') and (segname[i] == 'DNA1') and "
                          "(resid[i] <= 166) and (resid[i] >=  26))")  # NCP-1
        c11_filter.append("((name[i] == 'P') and (segname[i] == 'DNA1') and "
                          "(resid[i] <= 333) and (resid[i] >= 193))")  # NCP-2
        c11_filter.append("((name[i] == 'P') and (segname[i] == 'DNA1') and "
                          "(resid[i] <= 499) and (resid[i] >= 359))")  # NCP-3
        c11_filter.append("((name[i] == 'P') and (segname[i] == 'DNA1') and "
                          "(resid[i] <= 667) and (resid[i] >= 527))")  # NCP-4        
    
        all_ncps = []
    
        for i in xrange(len(c11_filter)):
            # load in the nucleosome with the linker protein
            linker = sasmol.SasMol(0)
            linker.read_pdb('gH5_NCP.pdb')
            # names = ['NCP1.pdb', 'NCP2.pdb', 'NCP3.pdb', 'NCP4.pdb']
            # for i in xrange(len(names)):
                # linker.read_pdb(names[i])
                # linker.setResname([name.replace('DA', 'ADE').replace(
                #'DC', 'CYT').replace('DG', 'GUA').replace('DT', 'THY')
                                   # for name in linker.resname()])
                # linker.write_pdb(names[i], 0, 'w')
                # print 'renamed resnames in %s' %names[i]
            # linker.write_pdb('gH5_NCP.pdb', 0, 'w')
        
            linker.center(0)
            # create the filter to select the coordinates for aligning
            # linker_filter = "((moltype[i] == 'dna' or moltype[i] == 'rna') and name[i] == 'P' and chain[i] == 'I' and resid[i] >= 14 and resid[i] <= 154)"
            linker_filter = "((moltype[i] == 'dna' or moltype[i] == 'rna') and name[i] == 'P' and chain[i] == 'I' and resid[i] >= 14 and resid[i] <= 64)"
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
    
            assert s_link == s_c11, "filter for c11 has %d atoms but H5 has %d atoms" % (s_c11, s_link)
            linker.align(0, coor_sub_2, com_sub_2, coor_sub_1, com_sub_1)
            # sub_1: the molecule to align too
            # sub_2: the molecule to be aligned
    
            # write the output to a new file
            ncp_filename = "NCP%d.pdb" % (i + 1)
            linker.write_pdb(ncp_filename, 0, 'w')
    
            all_ncps.append(linker)
    
        # ncp_array = combine_sasmols(all_ncps, combine_segnames=['I', 'J'])
        ncp_array = combine_sasmols(all_ncps)
        resname = [name.replace('DA', 'ADE').replace('DC', 'CYT').replace('DG', 'GUA').replace('DT', 'THY') for 
                   name in ncp_array.resname()]
        ncp_array.setResname(resname)
        ncp_array.write_pdb(pdb, 0, 'w') 
    
    bps = np.array([np.linspace(0, 167, 168), np.linspace(168, 1, 168)]).T
    ncp_dna_resids = [bps[[14, 154]]]*4
    dna_ids = [['1I', '1J'], ['2I', '2J'], ['3I', '3J'], ['4I', '4J']]
    ncp_dyad_resids = [bps[84], bps[84], bps[84], bps[84]]
    ncp_ref_atom_resids = [25, 25, 25, 25]    
    
    all_ncp_plot_vars, all_ncp_axes, all_ncp_origins = ta.get_tetramer_axes(
        pdb, ncp_dna_resids, dna_ids, ncp_dyad_resids, ncp_ref_atom_resids,
        array=ncp_array)
    phi, dxyz, plot_title = ta.get_tetramer_angles(all_ncp_axes, all_ncp_origins)
    geometry.show_ncps(all_ncp_plot_vars, title=plot_title)

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

def construct_ncp_array(ncp, phi, dxyz, dna_segnames, ncp_dna_resids, 
                        dyad_resids, ref_atom_resid, link_vars, pre_suf_vars,
                        save_name=None, adjust_ncp=None):
    '''
    given a list of sasmol objects, this will combine them into one 
    sasmol object
    
    inputs:
        ncp          : template ncp, either a sasmol object or filename for a pdb 
        phi          : orientation angle
        dxyz         : translation distances
        dna_segnames : ...
        
        save_name : optional input for what to save the result as 
        
    outputs:
        array     : resulting ncp array
        
    see also:
        align_gH5_to_c11
    '''
    assert len(phi) == len(dxyz), (
    'ERROR: inconsistent inputs, unclear definition parameters')
    n_ncp = len(phi)
    
    # load in NCP1 if it is not already a sasmol object
    if isinstance(ncp, basestring):
        filename = ncp
        ncp = sasmol.SasMol(0)
        ncp.read_pdb(filename)
    
    ncp_list = [ncp]
    
    # setup the masks for getting the axes and origin 
    # (should be consistent regardless coordinates)
    ncp_basis = ('(((segname[i] =="%s" and resid[i] >= %d and resid[i] <= %d ) or'
                 '  (segname[i] =="%s" and resid[i] <= %d and resid[i] >= %d) ) and'
                 '  name[i] == "C1\'")'
                 % (dna_segnames[0], ncp_dna_resids[0, 0], ncp_dna_resids[0, 1], 
                    dna_segnames[1], ncp_dna_resids[1, 0], ncp_dna_resids[1, 1]))
    error, ncp_mask = ncp.get_subset_mask(ncp_basis)
    
    dyad_basis = ('( segname[i] == "%s" and resid[i] == %d ) or'
                  '( segname[i] == "%s" and resid[i] == %d )'
                  % (dna_segnames[0], dyad_resids[0], 
                     dna_segnames[1], dyad_resids[1]) )
    error, dyad_mask = ncp.get_subset_mask(dyad_basis)
    
    ref_atom_basis = b2p.parse_basis("segname %s and resid %d and name C1\'"
                                     % (dna_segnames[0], ref_atom_resid))
    error, ref_atom_mask = ncp.get_subset_mask(ref_atom_basis)
    
    copy_mask = np.ones(ncp_mask.shape)
    ncp_origins = []
    ncp_axes = []
    for i in xrange(n_ncp):
        # get the axes and origin for NCP1
        ncp1 = ncp_list[i]
        if i == 0:
            debug=False
            (ncp1_origin, ncp1_axes, ncp_opt_params, ncp_dyad_mol, 
             ncp_plot_vars) = geometry.get_ncp_origin_and_axes(
                 ncp_mask, dyad_mask, dna_segnames, ncp1, ref_atom_mask,
                 debug=debug)
            ncp_origins.append(ncp1_origin)
            ncp_axes.append(ncp1_axes)
        else:
            ncp1_origin = ncp_origins[i]
            ncp1_axes = ncp_axes[i]
        
        # copy NCP1 to NCP2
        ncp2 = sasmol.SasMol(0)
        ncp_list.append(ncp2)
        error = ncp1.copy_molecule_using_mask(ncp2, copy_mask, 0)
        coor = ncp2.coor()[0]
        axes_coor = ncp1_axes + ncp1_origin
        all_coor = np.concatenate((ncp1_origin.reshape(1,3), axes_coor, coor))
        
        # rotate NCP2 about Z_1 by -phi_z
        all_coor = geometry.rotate_about_v(all_coor, ncp1_axes[2], -phi[i, 2])
        
        # rotate NCP2 about Y_1 by -phi_y
        all_coor = geometry.rotate_about_v(all_coor, ncp1_axes[1], -phi[i, 1])

        # rotate NCP2 about X_1 by -phi_x
        all_coor = geometry.rotate_about_v(all_coor, ncp1_axes[0], -phi[i, 0])
        
        # translate NCP2 a distance dx along X_1
        all_coor += dxyz[i, 0] * ncp1_axes[0]
        
        # translate NCP2 a distance dy along Y_1
        all_coor += dxyz[i, 1] * ncp1_axes[1]
    
        # translate NCP2 a distance dz along Z_1
        all_coor += dxyz[i, 2] * ncp1_axes[2]

        ncp_origins.append(all_coor[0])
        ncp_axes.append(all_coor[1:4] - all_coor[0])
        ncp2.setCoor(np.array([all_coor[4:]]))

    if adjust_ncp:
        reorient_ncps(ncp_list, ncp_axes, ncp_origins, adjust_ncp)

    array = combine_sasmols(ncp_list)

    if not save_name:
        save_name = 'tmp.pdb'
    array.write_pdb(save_name, 0, 'w')
    array.read_pdb(save_name)
    
    # align the linker dna
    ref_linker = sasmol.SasMol(0)
    ref_linker.read_pdb(link_vars.pdb)
    link_segnames = link_vars.segnames
    link_resids = link_vars.resids
    ncp_resids = link_vars.ncp_resids
    all_mols = [array]
    in_vars = inputs()
    in_vars.aa_goal = array
    linker_keep_basis = b2p.parse_basis('(segname DNA1 or segname DNA2) and '
                                        'resid >= %d and resid <= %d ' %
                                        (link_vars.keep[0], link_vars.keep[1]))
    error, keep_linker_mask = ref_linker.get_subset_mask(linker_keep_basis)
    for i in xrange(n_ncp):
        linker = sasmol.SasMol(0)
        error, mask = ref_linker.get_subset_mask('all')
        error = ref_linker.copy_molecule_using_mask(linker, mask, 0)
        in_vars.aa_move = linker
        move_vmd_basis = ('( segname %s and ( ( name C1\' and (resid %d or '
                          ' resid %d) ) or ( (name C2\' or name C1\') and '
                          'resid %d) ) )' % (link_segnames[0], link_resids[0,0],
                          link_resids[1,0], link_resids[2,0]) )
        in_vars.move_basis = b2p.parse_basis(move_vmd_basis)
    
        ncp_segnames = [[dna_segnames[0].replace('1', str(2*i+1)), 
                        dna_segnames[0].replace('1', str(2*i+2))],
                        [dna_segnames[0].replace('1', str(2*(i+1)+1)), 
                         dna_segnames[0].replace('1', str(2*(i+1)+2))]]  
        goal_vmd_basis = ('(segname %s and name C1\' and (resid %d or resid %d)'
                          ') or (segname %s and resid %d and (name C1\' or name'
                          ' C2\') )') % (ncp_segnames[0][0], ncp_resids[0,0], 
                          ncp_resids[1,0], ncp_segnames[1][0], ncp_resids[2,0])
        in_vars.goal_basis = b2p.parse_basis(goal_vmd_basis)

        align.align_mol(in_vars)
        keep_linker = sasmol.SasMol(0)
        error = in_vars.aa_move.copy_molecule_using_mask(
            keep_linker, keep_linker_mask, 0)
        new_segnames = [name.replace(link_segnames[0],'DNA0').replace(
            link_segnames[1],'DNA9') for name in keep_linker.segname()]
        keep_linker.setSegname(new_segnames)  
        all_mols.append(keep_linker)

    # align the prefix and suffix DNA
    segnames = pre_suf_vars.segnames

    # prefix #
    prefix = sasmol.SasMol(0)
    prefix.read_pdb(pre_suf_vars.pre_pdb)
    in_vars.aa_move = prefix
    pre_align_id = pre_suf_vars.pre_align_id
    move_vmd_basis = ('(segname %s and (name C1\' or name C2\') and '
                      '(resid %d or resid %d or resid %d) ) or '
                      '(segname %s and name C1\' and '
                      '(resid %d or resid %d or resid %d) )' % (
                      segnames[0], pre_align_id[0,0], 
                      pre_align_id[1,0], pre_align_id[2,0],
                      segnames[1], pre_align_id[0,1], 
                      pre_align_id[1,1], pre_align_id[2,1]))
    in_vars.move_basis = b2p.parse_basis(move_vmd_basis)
    
    pre_ref_id = pre_suf_vars.pre_ref_id
    goal_vmd_basis = ('(segname %s and (name C1\' or name C2\') and '
                      '(resid %d or resid %d or resid %d) ) or '
                      '(segname %s and name C1\' and '
                      '(resid %d or resid %d or resid %d) )' %
                      (dna_segnames[0], pre_ref_id[0,0], 
                      pre_ref_id[1,0], pre_ref_id[2,0],
                      dna_segnames[1], pre_ref_id[0,1], 
                      pre_ref_id[1,1], pre_ref_id[2,1]))
    in_vars.goal_basis = b2p.parse_basis(goal_vmd_basis)
    align.align_mol(in_vars)
    
    keep_pre = sasmol.SasMol(0)
    keep_pre_basis = ('resid <= %d or resid >= %d' % (pre_suf_vars.pre_keep[0],
                                                      pre_suf_vars.pre_keep[1]))
    error, keep_pre_mask = prefix.get_subset_mask(b2p.parse_basis(
        keep_pre_basis))
    error = in_vars.aa_move.copy_molecule_using_mask(keep_pre, keep_pre_mask, 0)
    new_segnames = [name.replace(segnames[0],'DNA0').replace(segnames[1],'DNA9') 
                    for name in keep_pre.segname()]
    keep_pre.setSegname(new_segnames)
    all_mols.append(keep_pre)
    
    # suffix #
    suffix = sasmol.SasMol(0)
    suffix.read_pdb(pre_suf_vars.suf_pdb)
    in_vars.aa_move = suffix
    suf_align_id = pre_suf_vars.suf_align_id
    move_vmd_basis = ('(segname %s and (name C1\' or name C2\') and '
                      '(resid %d or resid %d or resid %d) ) or '
                      '(segname %s and name C1\' and '
                      '(resid %d or resid %d or resid %d) )' % (
                      segnames[0], suf_align_id[0,0], 
                      suf_align_id[1,0], suf_align_id[2,0],
                      segnames[1], suf_align_id[0,1], 
                      suf_align_id[1,1], suf_align_id[2,1]))
    in_vars.move_basis = b2p.parse_basis(move_vmd_basis)
    
    suf_ref_id = pre_suf_vars.suf_ref_id
    goal_vmd_basis = ('(segname %s and (name C1\' or name C2\') and '
                      '(resid %d or resid %d or resid %d) ) or '
                      '(segname %s and name C1\' and '
                      '(resid %d or resid %d or resid %d) )' %
                      (ncp_segnames[1][0], suf_ref_id[0,0], 
                      suf_ref_id[1,0], suf_ref_id[2,0],
                      ncp_segnames[1][1], suf_ref_id[0,1], 
                      suf_ref_id[1,1], suf_ref_id[2,1]))
    in_vars.goal_basis = b2p.parse_basis(goal_vmd_basis)
    align.align_mol(in_vars)
    
    keep_suf = sasmol.SasMol(0)
    keep_suf_basis = ('resid >= %d or resid <= %d' % (pre_suf_vars.suf_keep[0],
                                                      pre_suf_vars.suf_keep[1]))
    error, keep_suf_mask = suffix.get_subset_mask(b2p.parse_basis(
        keep_suf_basis))
    error = in_vars.aa_move.copy_molecule_using_mask(keep_suf, keep_suf_mask, 0)
    new_segnames = [name.replace(segnames[0],'DNA0').replace(segnames[1],'DNA9') 
                    for name in keep_suf.segname()]
    keep_suf.setSegname(new_segnames)    
    all_mols.append(keep_suf)
    
    complete = combine_sasmols(all_mols)
    complete.write_pdb(save_name, 0, 'w')
    return complete

def reorient_ncps(ncp_list, ncp_axes, ncp_origins, adjust_ncp):
    method = adjust_ncp.method
    if 'spread' in method:
        # determine the rotation axes 
        v21 = ncp_origins[0] - ncp_origins[1]
        v23 = ncp_origins[2] - ncp_origins[1]
        v32 = -v23
        v34 = ncp_origins[3] - ncp_origins[2]
        ax2 = np.cross(v23, v21)
        ax3 = np.cross(v32, v34)
        spread_angle = adjust_ncp.spread_angle
        ncp2_origin = ncp_origins[1]
        ncp3_origin = ncp_origins[2]
        
        # rotate NCP1
        axes1_coor = ncp_axes[0] + ncp2_origin
        all_coor1 = np.concatenate((ncp2_origin.reshape(1,3), axes1_coor, 
                                    ncp_list[0].coor()[0]))
        all_coor1 = geometry.rotate_about_v(all_coor1, ax2, spread_angle)

        # rotate NCP2
        axes2_coor = ncp_axes[1] + ncp2_origin
        all_coor2 = np.concatenate((ncp2_origin.reshape(1,3), axes2_coor, 
                                    ncp_list[1].coor()[0]))
        all_coor2 = geometry.rotate_about_v(all_coor2, ax2, spread_angle/2.0)
        
        # store the output
        assert 6 == np.isclose((ncp2_origin, ncp2_origin), 
                          (all_coor1[0], all_coor1[0])).sum(), (
                              'ERROR: origin shifted') # sum(True) = 1
        ncp_axes[0] = all_coor1[1:4] - all_coor1[0]
        ncp_list[0].setCoor(np.array([all_coor1[4:]]))
        ncp_axes[1] = all_coor2[1:4] - all_coor2[0]
        ncp_list[1].setCoor(np.array([all_coor2[4:]]))

        # rotate NCP3
        axes3_coor = ncp_axes[2] + ncp3_origin
        all_coor3 = np.concatenate((ncp3_origin.reshape(1,3), axes3_coor, 
                                    ncp_list[2].coor()[0]))
        all_coor3 = geometry.rotate_about_v(all_coor3, ax3, spread_angle)
    
        # rotate NCP4
        axes4_coor = ncp_axes[3] + ncp3_origin
        all_coor4 = np.concatenate((ncp3_origin.reshape(1,3), axes4_coor, 
                                    ncp_list[3].coor()[0]))
        all_coor4 = geometry.rotate_about_v(all_coor4, ax3, spread_angle)
    
        # store the output
        assert 6 == np.isclose((ncp3_origin, ncp3_origin), 
                               (all_coor3[0], all_coor3[0])).sum(), (
                                   'ERROR: origin shifted') # sum(True) = 1
        ncp_axes[2] = all_coor3[1:4] - all_coor3[0]
        ncp_list[2].setCoor(np.array([all_coor3[4:]]))
        ncp_axes[3] = all_coor4[1:4] - all_coor4[0]
        ncp_list[3].setCoor(np.array([all_coor4[4:]]))
    
    if 'slide' in method:
        mag = adjust_ncp.mag
        for i in xrange(len(adjust_ncp.mv_ncp)):
            mv_ncp = adjust_ncp.mv_ncp[i]
            i_axes = adjust_ncp.i_axes[i]
            ax1 = ncp_axes[i_axes[0][0]][i_axes[0][1]]
            ax2 = ncp_axes[i_axes[1][0]][i_axes[1][1]] 
            stack_axis = (ax1 + ax2)
            stack_axis = stack_axis/np.sqrt(stack_axis.dot(stack_axis)) # unit v
            displacement = mag * stack_axis
            ncp_list[mv_ncp[0]].setCoor(ncp_list[mv_ncp[0]].coor() - displacement)
            ncp_list[mv_ncp[1]].setCoor(ncp_list[mv_ncp[1]].coor() + displacement)

    if 'rotate_sym' in method:
        center = np.array(ncp_origins).mean(axis=0)
        angles = [adjust_ncp.angle[0], -adjust_ncp.angle[1]]
        for (i, mv_ncp) in enumerate(adjust_ncp.mv_ncp):
            mv_ncp = adjust_ncp.mv_ncp[i] 
            # # ~~~ the axis halfway between the two cylinder axes ~~~ #
            # i_axes = adjust_ncp.i_axes[i]
            # ax1 = ncp_axes[i_axes[0][0]][i_axes[0][1]]
            # ax2 = ncp_axes[i_axes[1][0]][i_axes[1][1]] 
            # slide_axis = (ax1 + ax2)
            # # ~~~ the axis passing through the center of both NCPs ~~~ # 
            stack_axis = ncp_origins[mv_ncp[1]] - ncp_origins[mv_ncp[0]]

            stack_axis = stack_axis/np.sqrt(stack_axis.dot(stack_axis)) # unit v
            for i in xrange(2):
                ncp = ncp_list[mv_ncp[i]]
                axes_coor = ncp_axes[mv_ncp[i]] + center
                coor = ncp.coor()[0]
                all_coor = np.concatenate((center.reshape(1,3), axes_coor, coor))
                
                # perform the rotation
                all_coor = geometry.rotate_about_v(all_coor, stack_axis, angles[i])
    
                # store the output
                assert all(np.isclose(center, all_coor[0])), (
                    'ERROR: origin shifted, check input')
                ncp_axes[mv_ncp[i]] = all_coor[1:4] - all_coor[0]
                ncp.setCoor(np.array([all_coor[4:]]))            
    
    if 'rotate_asym' in method:
        center = np.array(ncp_origins).mean(axis=0)
        for (i, mv_ncp) in enumerate(adjust_ncp.mv_ncp):
            angles = [adjust_ncp.angle[i], -adjust_ncp.angle[i]]
            # # ~~~ the axis halfway between the two cylinder axes ~~~ #
            # i_axes = adjust_ncp.i_axes[i]
            # ax1 = ncp_axes[i_axes[0][0]][i_axes[0][1]]
            # ax2 = ncp_axes[i_axes[1][0]][i_axes[1][1]] 
            # slide_axis = (ax1 + ax2)
            # # ~~~ the axis passing through the center of both NCPs ~~~ # 
            stack_axis = ncp_origins[mv_ncp[1]] - ncp_origins[mv_ncp[0]]

            stack_axis = stack_axis/np.sqrt(stack_axis.dot(stack_axis)) # unit v
            for (j, i_ncp) in enumerate(mv_ncp):
                ncp = ncp_list[i_ncp]
                axes_coor = ncp_axes[i_ncp] + center
                coor = ncp.coor()[0]
                all_coor = np.concatenate((center.reshape(1,3), axes_coor, coor))
                
                # perform the rotation
                all_coor = geometry.rotate_about_v(all_coor, stack_axis, angles[j])
    
                # store the output
                assert all(np.isclose(center, all_coor[0])), (
                    'ERROR: origin shifted, check input')
                ncp_axes[i_ncp] = all_coor[1:4] - all_coor[0]
                ncp.setCoor(np.array([all_coor[4:]]))            
        
    if 'twist' in method:
        mv_ncp = adjust_ncp.mv_ncp
        i_axes = adjust_ncp.i_axes # [which ncp, which axis]
        angles = adjust_ncp.angles
        origin = adjust_ncp.origin
        
        assert len(i_axes) == len(angles) == len(mv_ncp) == len(origin), (
            'ERROR: mismatch in input')
        
        for (i, i_ncp) in enumerate(mv_ncp):
            
            # setup the rotation
            ncp = ncp_list[i_ncp]
            ax2 = ncp_axes[i_axes[i][0]][i_axes[i][1]]
            
            rotation_origin = ncp_origins[origin[i]]
            axes_coor = ncp_axes[i_ncp] + rotation_origin
            coor = ncp.coor()[0]
            all_coor = np.concatenate((rotation_origin.reshape(1,3), axes_coor, 
                                       coor))
    
            # perform the rotation
            all_coor = geometry.rotate_about_v(all_coor, ax2, angles[i])
            
            # store the output
            assert all(np.isclose(rotation_origin, all_coor[0])), (
                'ERROR: origin shifted')
            ncp_axes[i_ncp] = all_coor[1:4] - all_coor[0]
            ncp.setCoor(np.array([all_coor[4:]]))

if __name__ == '__main__':
    # align_gH5_to_c11()
    
    # ncp = 'NCP1.pdb'
    # dna_segnames = ['I', 'J']
    # w601 = [14, 154]
    # bps = np.array([np.linspace(0, 167, 168), np.linspace(168, 1, 168)]).T
    ncp = 'gH5_1x164.pdb'
    dna_segnames = ['DNA1', 'DNA2']
    w601 = [12, 152]
    ncp_link_match = [163, 164, 1, 2]
    bps = np.array([np.linspace(0, 164, 165), np.linspace(165, 1, 165)]).T
    ref_atom_resid = 23
    
    link_vars = inputs()
    link_vars.pdb = 'linker.pdb'
    link_vars.segnames = ['DNA1', 'DNA2']
    link_vars.resids = np.array([[1, 7], [2, 6], [7, 1]])
    link_vars.keep = [3,5]
    link_vars.ncp_resids = bps[ncp_link_match]
    
    pre_suf_vars = inputs()
    pre_suf_vars.segnames = ['DNA1', 'DNA2']

    # pre_suf_vars.pre_pdb = 'dna_prefix190.pdb'
    # save_name = 'complete_gH5x4_190.pdb'
    # pre_suf_vars.pre_pdb = 'dna_prefix191.pdb'
    # save_name = 'complete_gH5x4_191.pdb'
    pre_suf_vars.pre_pdb = 'dna_prefix218.pdb'
    save_name = 'complete_gH5x4_218.pdb'
    pre_bps = np.array([np.linspace(0, 17, 18), np.linspace(694, 677, 18)]).T
    pre_suf_vars.pre_keep = pre_bps[14]
    pre_suf_vars.pre_align_id = pre_bps[-3:]
    pre_suf_vars.pre_ref_id = bps[1:4]


    pre_suf_vars.suf_pdb = 'dna_suffix.pdb'
    suf_bps = np.array([np.linspace(694, 677, 18), np.linspace(0, 17, 18)]).T
    pre_suf_vars.suf_keep = suf_bps[14]
    pre_suf_vars.suf_align_id = suf_bps[-3:]
    pre_suf_vars.suf_ref_id = bps[-3:]
    
    phi_file = 'gH5c11_r_phi.txt'
    dxyz_file = 'gH5c11_r_dxyz.txt'
    phi = np.loadtxt(phi_file)
    dxyz = np.loadtxt(dxyz_file)

    adjust_ncp = None
    # adjust_ncp = inputs()
    # adjust_ncp.mv_ncp = [1,3]
    # adjust_ncp.i_axes = [[1,0],[3,0]]
    # # adjust_ncp.i_axes = [[1,0],[1,0]]
    # adjust_ncp.angles = [20, 20]
    save_name = '150602_gH5x4.pdb'

    ncp_dna_resids = bps[[w601[0], w601[1]]]
    dyad_resids = bps[(w601[1] - w601[0])/2 + w601[0]]
    array = construct_ncp_array(ncp, phi, dxyz, dna_segnames, ncp_dna_resids, 
                                dyad_resids, ref_atom_resid, link_vars, 
                                pre_suf_vars, save_name = save_name, 
                                adjust_ncp = adjust_ncp)
    # array.write_pdb('complete_gH5x4.pdb', 0, 'w')
    print '\m/ >.< \m/'