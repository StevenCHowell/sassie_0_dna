#!/usr/bin/env python
# coding:utf-8
"""
  Author:  Steven C. Howell
  Purpose: construct gH5 nucleosome arrays incrementally twisting NCP 2 and 4
  Created: 5 May 2015
  $Id: build_4x167_gH5.py 208 2015-05-04 18:45:41Z xraylab $
00000000011111111112222222222333333333344444444445555555555666666666677777777778
12345678901234567890123456789012345678901234567890123456789012345678901234567890
"""

import numpy as np
import x_dna.build_mol.array_from_angles.build_4x167_gH5 as build

ncp = 'gH5_1x164.pdb'
dna_segnames = ['DNA1', 'DNA2']
w601 = [12, 152]
ncp_link_match = [163, 164, 1, 2]
bps = np.array([np.linspace(0, 164, 165), np.linspace(165, 1, 165)]).T
ref_atom_resid = 23

link_vars = build.inputs()
link_vars.pdb = 'linker.pdb'
link_vars.segnames = ['DNA1', 'DNA2']
link_vars.resids = np.array([[1, 7], [2, 6], [7, 1]])
link_vars.keep = [3, 5]
link_vars.ncp_resids = bps[ncp_link_match]

pre_suf_vars = build.inputs()
pre_suf_vars.segnames = ['DNA1', 'DNA2']

pre_suf_vars.pre_pdb = 'dna_prefix218.pdb'
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

adjust_ncp = build.inputs()
adjust_ncp.mv_ncp = [1, 3]
adjust_ncp.origin = [1, 3]
adjust_ncp.method = 'twist'
ncp_dna_resids = bps[[w601[0], w601[1]]]
dyad_resids = bps[(w601[1] - w601[0]) / 2 + w601[0]]

for i in xrange(-1, 4):
    adjust_ncp.angles = [i, i]
    adjust_ncp.i_axes = [[1, 0], [3, 0]]
    save_name = 'gH5x4_%s_n2n4_r24.pdb' % str(i).replace('-', 'm')
    array = build.construct_ncp_array(ncp, phi, dxyz, dna_segnames,
                                      ncp_dna_resids, dyad_resids,
                                      ref_atom_resid, link_vars,
                                      pre_suf_vars, save_name=save_name,
                                      adjust_ncp=adjust_ncp)

    # adjust_ncp.i_axes = [[1,0],[1,0]]
    # save_name = 'gH5x4_%s_n2n4_r22.pdb' % str(i).replace('-','m')
    # array = build.construct_ncp_array(ncp, phi, dxyz, dna_segnames,
    # ncp_dna_resids, dyad_resids,
    # ref_atom_resid, link_vars,
    # pre_suf_vars, save_name = save_name,
    # adjust_ncp = adjust_ncp)

print '\m/ >.< \m/'
