#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose:
# Created: 
#
# $Id: $
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

# import ncp_angles as na
import sassie.sasmol.sasmol as sasmol
import numpy as np
import sassie_1_na.util.geometry as geometry
import sassie_1_na.util.basis_to_python as basis_to_python
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.core.umath_tests import inner1d
try: 
    import cPickle as pickle
except:
    import pickle as pickle
    
def get_tetramer_axes(pdb, ncp_dna_resids, dna_ids, ncp_dyad_resids):
    pkl_file = pdb[:-3] + 'pkl'
    try:
        pkl_in = open(pkl_file, 'rb')
        load_masks = True
        all_ncp_bases      = pickle.load(pkl_in)
        all_ncp_masks      = pickle.load(pkl_in)
        all_dyad_bases     = pickle.load(pkl_in)
        all_dyad_masks     = pickle.load(pkl_in)
        all_ncp_origins    = pickle.load(pkl_in)
        all_ncp_axes       = pickle.load(pkl_in)
        all_ncp_opt_params = pickle.load(pkl_in)
        all_ncp_dyad_mol   = pickle.load(pkl_in)
        all_ncp_plot_vars  = pickle.load(pkl_in)
        pkl_in.close()
        print 'using masks from: %s' % pkl_file
    except:
        if 'load_masks' in locals():
            pkl_in.close()
        load_masks = False
        n_ncps = len(ncp_dna_resids)
        all_ncp_bases      = [None] * n_ncps
        all_ncp_masks      = [None] * n_ncps
        all_dyad_bases     = [None] * n_ncps
        all_dyad_masks     = [None] * n_ncps
        all_ncp_origins    = [None] * n_ncps
        all_ncp_axes       = [None] * n_ncps
        all_ncp_opt_params = [None] * n_ncps
        all_ncp_dyad_mol   = [None] * n_ncps 
        all_ncp_plot_vars  = [None] * n_ncps
        print 'creating masks from input variables'

    array = sasmol.SasMol(0)
    array.read_pdb(pdb)

    # # re-orient the array
    # coor = array.coor()
    # coor[0] = geometry.transform_coor(coor[0], np.array([1, 0, 0]), np.array([0, 0, 0]))
    # array.setCoor(coor)
    # array.write_pdb('gH5c11_r.pdb', 0, 'w')

    errors = []
    for (i, resids) in enumerate(ncp_dna_resids):
        print 'fitting NCP %d' % (i+1)
        if load_masks:
            ncp_mask       = all_ncp_masks[i]
            dyad_mask      = all_dyad_masks[i]
            ncp_opt_params = all_ncp_opt_params[i]
        else:
            ncp_basis_vmd = ("((segname %s and resid >= %d and resid <=  %d) or"
                             " (segname %s and resid <= %d and resid >= %d) ) and name C1' " 
                             % (dna_ids[i][0], resids[0,0], resids[1,0], dna_ids[i][1], resids[0,1], resids[1,1]))
            ncp_basis = basis_to_python.parse_basis(ncp_basis_vmd)
            error, ncp_mask = array.get_subset_mask(ncp_basis)
            errors.append(errors)
            
            dyad_basis = ('( segname[i] == "%s" and resid[i] == %d ) or'
                          '( segname[i] == "%s" and resid[i] == %d )'
                          % (dna_ids[i][0], ncp_dyad_resids[i][0], dna_ids[i][1], ncp_dyad_resids[i][1]) )
            error, dyad_mask = array.get_subset_mask(dyad_basis)
            errors.append(errors)
            ncp_opt_params = None
            all_ncp_bases[i]      = ncp_basis
            all_ncp_masks[i]      = ncp_mask
            all_dyad_bases[i]     = dyad_basis
            all_dyad_masks[i]     = dyad_mask

        ncp_origin, ncp_axes, ncp_opt_params, ncp_dyad_mol, ncp_plot_vars = geometry.get_ncp_origin_and_axes(
            ncp_mask, dyad_mask, dna_ids[i], array, prev_opt_params=ncp_opt_params, debug=True)
        
        # store all the variables for debug purposes
        all_ncp_origins[i]    = ncp_origin
        all_ncp_axes[i]       = ncp_axes
        all_ncp_opt_params[i] = ncp_opt_params
        all_ncp_dyad_mol[i]   = ncp_dyad_mol
        all_ncp_plot_vars[i]  = ncp_plot_vars

    print 'saving masks and other results to: %s' % pkl_file
    pkl_out = open(pkl_file, 'wb')
    pickle.dump(all_ncp_bases, pkl_out, -1)
    pickle.dump(all_ncp_masks, pkl_out, -1)
    pickle.dump(all_dyad_bases, pkl_out, -1)
    pickle.dump(all_dyad_masks, pkl_out, -1)
    pickle.dump(all_ncp_origins, pkl_out, -1)
    pickle.dump(all_ncp_axes, pkl_out, -1)
    pickle.dump(all_ncp_opt_params, pkl_out, -1)
    pickle.dump(all_ncp_dyad_mol, pkl_out, -1)
    pickle.dump(all_ncp_plot_vars, pkl_out, -1)
    pkl_out.close()
    return all_ncp_plot_vars, all_ncp_axes, all_ncp_origins
    
def main():
    NotImplemented


if __name__ == '__main__':
    pdb = 'gH5c11_r.pdb'  # rotated to be more perpendicular to the x-y plane
    # pdb = 'gH5c11.pdb' # this creates a problem where the z-axes are anti-aligned
    bps = np.array([np.linspace(0,694,695), np.linspace(694, 0, 695)]).T

    dna_ids = [['DNA1', 'DNA2']]*4
    ncp_dna_resids = [bps[[26, 166]], bps[[193, 333]], bps[[360, 500]], bps[[527, 667]]]
    ncp_dyad_resids = [bps[96], bps[263], bps[430], bps[597]]

    # get tetramer axes
    pkl_file = pdb[:-3] + 'pkl'
    try:
        pkl_in = open(pkl_file, 'rb')
        load_masks = True
        all_ncp_bases      = pickle.load(pkl_in)
        all_ncp_masks      = pickle.load(pkl_in)
        all_dyad_bases     = pickle.load(pkl_in)
        all_dyad_masks     = pickle.load(pkl_in)
        all_ncp_origins    = pickle.load(pkl_in)
        all_ncp_axes       = pickle.load(pkl_in)
        all_ncp_opt_params = pickle.load(pkl_in)
        all_ncp_dyad_mol   = pickle.load(pkl_in)
        all_ncp_plot_vars  = pickle.load(pkl_in)
        pkl_in.close()
    except:
        all_ncp_plot_vars, all_ncp_axes, all_ncp_origins = get_tetramer_axes(pdb, ncp_dna_resids, dna_ids, ncp_dyad_resids)
    geometry.show_ncps(all_ncp_plot_vars)
    
    if False:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        n1_origin = all_ncp_plot_vars[0].ncp_origin
        styles = ['r-','g-','b-']
        labels = ['NCP1_X', 'NCP1_Y', 'NCP1_Z']
        for (j, axis) in enumerate(all_ncp_plot_vars[0].ncp_axes):
            axes_vec = np.vstack((n1_origin, axis*15 + n1_origin))
            ax.plot(axes_vec[:,0], axes_vec[:,1], axes_vec[:,2], styles[j], label=labels[j])
        n2_x_axis = all_ncp_plot_vars[1].ncp_axes[0]
        n2_x_axis_vec = np.vstack((n1_origin, n2_x_axis*15 + n1_origin))
        ax.plot(n2_x_axis_vec[:,0], n2_x_axis_vec[:,1], n2_x_axis_vec[:,2], 'm-', label='NCP2_X')
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.auto_scale_xyz([0, 100], [0, 100], [0, 100])
        # plt.axis('equal')
        plt.legend(loc='upper left', numpoints=1, bbox_to_anchor=(1, 0.5))
        plt.show()
        
    all_ncp_axes = np.array(all_ncp_axes)
    all_ncp_origins = np.array(all_ncp_origins)
    
    # get tetramer angles
    n_ncp = len(ncp_dna_resids)
    all_ncp1_axes = all_ncp_axes[:-1] # ncp 1-3
    all_ncp2_axes = all_ncp_axes[ 1:] # ncp 2-4
    all_ncp1_origins = all_ncp_origins[:-1] # ncp 1-3
    all_ncp2_origins = all_ncp_origins[ 1:] # ncp 2-4

    # bending
    phi_r = np.arccos(inner1d(all_ncp1_axes[:, 0, :], all_ncp2_axes[:, 0, :]))
    phi_switch = inner1d(np.cross(all_ncp1_axes[:, 0, :], all_ncp2_axes[:, 0, :]), all_ncp1_axes[:, 2, :]) < 0
    phi_r[phi_switch] = 2 * np.pi - phi_r[phi_switch]
    phi_d = phi_r * 180 / np.pi
    
    # twist (Victor's group from NIH uses rho, my psi/2)
    psi_r = np.arccos(inner1d(all_ncp1_axes[:, 2, :], all_ncp2_axes[:, 2, :]))
    psi_switch = inner1d(np.cross(all_ncp2_axes[:, 2, :], all_ncp1_axes[:, 2, :]), all_ncp1_axes[:, 0, :]) < 0
    psi_r[psi_switch] = 2 * np.pi - psi_r[psi_switch]
    psi_d = psi_r * 180 / np.pi
    
    # rise     
    h = inner1d(all_ncp1_axes[:, 2, :], all_ncp2_origins - all_ncp1_origins)
    
    data = 'N4merH5TE_zeroCon.iq'
    # data = 'N4merH5TE_zeroCon.i0q'
    # data = 'N4merH5Mg1_zeroCon.iq'
    # data = 'N4merH5Mg1_zeroCon.i0q'

    plot_title = (r'$\phi$ (bend): (%0.1f, %0.1f, %0.1f); '
                  r'$\psi$ (twist): (%0.1f, %0.1f, %0.1f); '
                  r'h (rise): (%0.1f, %0.1f, %0.1f)'
                  % (phi_d[0], phi_d[1], phi_d[2], 
                     psi_d[0], psi_d[1], psi_d[2],
                     h[0], h[1], h[2]))
    geometry.show_ncps(all_ncp_plot_vars, title=plot_title)
    print '\m/ >.< \m/'    
    