#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Calculate the angles between the nucleosomes in a tetramer array
# Created: 15 April 2015
#
# $Id$
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#

# import ncp_angles as na
import sassie.sasmol.sasmol as sasmol
import numpy as np
import x_dna.util.geometry as geometry
import x_dna.util.basis_to_python as basis_to_python
# import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.core.umath_tests import inner1d
try: 
    import cPickle as pickle
except:
    import pickle as pickle
    
def get_tetramer_axes(pdb, ncp_dna_resids, dna_ids, ncp_dyad_resids, 
                      ncp_ref_atom_resids, array=None):
    pkl_file = pdb[:-3] + 'pkl'
    try:
        pkl_in = open(pkl_file, 'rb')
        load_masks = True
        all_ncp_bases      = pickle.load(pkl_in)
        all_ncp_masks      = pickle.load(pkl_in)
        all_dyad_bases     = pickle.load(pkl_in)
        all_dyad_masks     = pickle.load(pkl_in)
        all_ref_masks      = pickle.load(pkl_in)
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
        all_ref_atom_masks = [None] * n_ncps
        all_ncp_origins    = [None] * n_ncps
        all_ncp_axes       = [None] * n_ncps
        all_ncp_opt_params = [None] * n_ncps
        all_ncp_dyad_mol   = [None] * n_ncps 
        all_ncp_plot_vars  = [None] * n_ncps
        print 'creating masks from input variables'

    # check if pdb is a filename or a sasmol object
    if not array:
        array = sasmol.SasMol(0)
        array.read_pdb(pdb)

    # # re-orient the array
    # coor = array.coor()
    # coor[0] = geometry.transform_coor(coor[0], np.array([1, 0, 0]), 
                                      # np.array([0, 0, 0]))
    # array.setCoor(coor)
    # array.write_pdb('gH5c11_r.pdb', 0, 'w')

    errors = []
    for (i, resids) in enumerate(ncp_dna_resids):
        print 'fitting NCP %d' % (i+1)
        if load_masks:
            ncp_mask       = all_ncp_masks[i]
            dyad_mask      = all_dyad_masks[i]
            ncp_opt_params = all_ncp_opt_params[i]
            ref_atom_mask  = all_ref_atom_masks[i]
        else:
            ncp_basis_vmd = ("((segname %s and resid >= %d and resid <=  %d) or"
                             " (segname %s and resid <= %d and resid >= %d) ) "
                             "and name C1' " 
                             % (dna_ids[i][0], resids[0,0], resids[1,0], 
                                dna_ids[i][1], resids[0,1], resids[1,1]))
            ncp_basis = basis_to_python.parse_basis(ncp_basis_vmd)
            error, ncp_mask = array.get_subset_mask(ncp_basis)
            errors.append(errors)
            
            dyad_basis = ('( segname[i] == "%s" and resid[i] == %d ) or'
                          '( segname[i] == "%s" and resid[i] == %d )'
                          % (dna_ids[i][0], ncp_dyad_resids[i][0], 
                             dna_ids[i][1], ncp_dyad_resids[i][1]) )
            error, dyad_mask = array.get_subset_mask(dyad_basis)
            errors.append(errors)

            ref_atom_basis_vmd = ("(segname %s and resid %d and name C1\')" %
                                  (dna_ids[i][0], ncp_ref_atom_resids[i]))
            ref_atom_basis = basis_to_python.parse_basis(ref_atom_basis_vmd)
            error, ref_atom_mask = array.get_subset_mask(ref_atom_basis)
            
            ncp_opt_params = None
            all_ncp_bases[i]      = ncp_basis
            all_ncp_masks[i]      = ncp_mask
            all_dyad_bases[i]     = dyad_basis
            all_dyad_masks[i]     = dyad_mask
            all_ref_atom_masks[i] = ref_atom_mask
            
        (ncp_origin, ncp_axes, ncp_opt_params, ncp_dyad_mol, ncp_plot_vars
         ) = geometry.get_ncp_origin_and_axes(ncp_mask, dyad_mask, dna_ids[i], 
                                              array, ref_atom_mask, debug=True,
                                              prev_opt_params=ncp_opt_params)
        
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
    pickle.dump(all_ref_atom_masks, pkl_out, -1)
    pickle.dump(all_ncp_origins, pkl_out, -1)
    pickle.dump(all_ncp_axes, pkl_out, -1)
    pickle.dump(all_ncp_opt_params, pkl_out, -1)
    pickle.dump(all_ncp_dyad_mol, pkl_out, -1)
    pickle.dump(all_ncp_plot_vars, pkl_out, -1)
    pkl_out.close()
    return all_ncp_plot_vars, all_ncp_axes, all_ncp_origins

def get_tetramer_angles(all_ncp_axes, all_ncp_origins):
    # get tetramer angles
    all_ncp_axes = np.array(all_ncp_axes)
    all_ncp_origins = np.array(all_ncp_origins)
    
    n_ncp = len(all_ncp_axes)
    ncp1_axes = all_ncp_axes[:-1] # ncp 1-3
    ncp2_axes = all_ncp_axes[ 1:] # ncp 2-4
    ncp1_origins = all_ncp_origins[:-1] # ncp 1-3
    ncp2_origins = all_ncp_origins[ 1:] # ncp 2-4

    #~ these angles and translations are different from Victor Zhurkin @ NIH ~#
    # phi_x: angle that aligns Z_2 to the X_1-Z_1 plane
    phi = np.zeros((n_ncp-1, 3))
    for i in xrange(n_ncp-1):
        phi[i] = geometry.get_alignment_angles(ncp1_axes[i], ncp2_axes[i])

    # R_align = geometry.rotate_v2_to_v1(all_ncp1_axes[:, 2, :], all_ncp2_axes[:, 2, :])
    # # np.einsum('ijk,ik->ij', R_align, all_ncp2_axes[:, 2, :])
    # # np.einsum('jk,lk->lj', R_align[0], all_ncp2_axes[0])
    # rot_ncp2_axes = np.einsum('ijk,ilk->ilj', R_align, all_ncp2_axes)
    
    # # bending
    # phi_r = np.arccos(inner1d(all_ncp1_axes[:, 0, :], rot_ncp2_axes[:, 0, :]))
    # phi_switch = inner1d(np.cross(all_ncp1_axes[:, 0, :], rot_ncp2_axes[:, 0, :]), all_ncp1_axes[:, 2, :]) < 0
    # phi_r[phi_switch] = 2 * np.pi - phi_r[phi_switch]
    # phi_d = phi_r * 180 / np.pi
    

    # # rise (not effective without the radius values)
    # h0 = inner1d(all_ncp1_axes[:, 2, :], all_ncp2_origins - all_ncp1_origins)
    
    # radius and rise
    # # solve this linear system: r1 * X_1 - r2 * X_2 + h * Z_1 = O_2 - O_1 (not effective)
    # results = []
    # for i in xrange(len(all_ncp1_origins)):
        # a = np.array([all_ncp1_axes[i,0,:], -all_ncp2_axes[i,0,:], all_ncp1_axes[i,2,:]]).transpose()
        # b = all_ncp2_origins[i] - all_ncp1_origins[i]
        # results.append(np.linalg.solve(a, b))
    # [r1, r2, h] = np.array(results).transpose()
    
    # solve this linear system: x * X_1 + y * Y_1 + z * Z_1 = O_2 - O_1 (reversible)
    results = []
    for i in xrange(len(ncp1_origins)):
        a = np.array([ncp1_axes[i,0,:], ncp1_axes[i,1,:], ncp1_axes[i,2,:]]).transpose()
        b = ncp2_origins[i] - ncp1_origins[i]
        results.append(np.linalg.solve(a, b))
    dxyz = np.array(results)
    
    
    for i in xrange(3):
        p1 = (ncp1_origins[i] + dxyz[i, 0] * ncp1_axes[i, 0, :] + 
              dxyz[i, 1] * ncp1_axes[i, 1, :] + dxyz[i, 2] * ncp1_axes[i, 2, :])
        p2 = ncp2_origins[i]
        assert np.allclose(p1, p2), 'WARNING: problem with results'

    plot_title = (r'($\phi_x$, $\phi_y$ $\phi_z)$: (%0.1f, %0.1f, %0.1f), (%0.1f, %0.1f, %0.1f), (%0.1f, %0.1f, %0.1f)'
                  '\n'
                  r'($d_x$, $d_y$, $d_z$): (%0.1f, %0.1f, %0.1f), (%0.1f, %0.1f, %0.1f), (%0.1f, %0.1f, %0.1f)'
                  % (phi[0,0],  phi[0,1],  phi[0,2], 
                     phi[1,0],  phi[1,1],  phi[1,2], 
                     phi[2,0],  phi[2,1],  phi[2,2], 
                     dxyz[0,0], dxyz[0,1], dxyz[0,2], 
                     dxyz[1,0], dxyz[1,1], dxyz[1,2], 
                     dxyz[2,0], dxyz[2,1], dxyz[2,2]))
    
    return phi, dxyz, plot_title
    
def main():
    NotImplemented


if __name__ == '__main__':
    bps = np.array([np.linspace(0,694,695), np.linspace(694, 0, 695)]).T
    dna_ids = [['DNA1', 'DNA2']]*4
    ncp_dna_resids = [bps[[26, 166]], bps[[193, 333]], bps[[360, 500]], bps[[527, 667]]]
    ncp_dyad_resids = [bps[96], bps[263], bps[430], bps[597]]
    ncp_ref_atom_resids = [37, 204, 371, 538]
    # the 1ZBB structure with corrected DNA sequence and completed proteins (default)
    pdb = 'gH5c11_r.pdb'  # rotated to be more perpendicular to the x-y plane
    ## angle between sequenctial NCP z-axes:  [ 93.7409857   93.38903002  93.79568914]
    ## angle between stacked NCP z-axes:      [ 12.58114363  13.1607299 ]
    ## opening angles between NCPs:           [ 22.8974572   23.02477343]
    ## stack distances:                       [57.807627493800133, 58.212257694966162]

    # the unrotated 1ZBB structure # this creates a problem where the z-axes are anti-aligned
    # pdb = '/home/schowell/data/code/pylib/x_dna/build_mol/array_from_angles/gH5c11.pdb'
    ## angle between sequenctial NCP z-axes:  [ 93.74106829  93.38933082  93.79580714]
    ## angle between stacked NCP z-axes:      [ 12.5809731   13.15984983]
    ## opening angles between NCPs:           [ 22.89739159  23.02455461]
    ## stack distances:                       [57.807471776623998, 58.211688477840951]
    
    # parameters for rigid body models
    bps = np.array([np.linspace(0,164,165), np.linspace(164, 0, 165)]).T
    dna_ids = [['DNA1', 'DNA2'],['DNA3', 'DNA4'],['DNA5', 'DNA6'],['DNA7', 'DNA8']]
    ncp_dna_resids = [bps[[12, 152]]]*4
    ncp_dyad_resids = [bps[82]]*4
    ncp_ref_atom_resids = [23]*4

    # the 2 degree opening adjustment using rigid body modification 
    # pdb = '/home/schowell/data/myData/manualStructures/gH5_opening/pdb/gH5x4_opening_2d.pdb'
    ## angle between sequenctial NCP z-axes:  [ 93.8890999   93.79122856  93.97284097]
    ## angle between stacked NCP z-axes:      [ 14.99871346  15.36996521]
    ## opening angles between NCPs:           [ 24.91176245  25.01137675]
    ## stack distances:                       [62.745230145841695, 63.144230926884156]
    
    # the -10 degree opening angle
    pdb = '/home/schowell/data/myData/manualStructures/gH5_opening/pdb/gH5x4_opening_m10d.pdb'
    ## angle between sequenctial NCP z-axes:  [ 94.52648864  92.37749212  93.97302923]
    ## angle between stacked NCP z-axes:      [ 2.29559141  2.83792858]
    ## opening angles between NCPs:           [ 12.90762848  13.01391677]
    ## stack distances:                       [32.972430152964371, 33.373636182110772]

    pdb = '/home/schowell/data/myData/manualStructures/gH5_opening/pdb/gH5x4_opening_10d.pdb'

    # the opening starting structure
    # pdb = '/home/schowell/data/myData/manualStructures/gH5_opening/pdb/gH5x4_opening_0d.pdb'
    ## angle between sequenctial NCP z-axes:  [ 94.00369639  93.60561237  93.97278084]
    ## angle between stacked NCP z-axes:      [ 12.62481159  13.19528317]
    ## opening angles between NCPs:           [ 22.91110313  23.01163443]
    ## stack distances:                       [57.809280925176395, 58.205898862508221]
    
    # # the -14 degree twist adjustment using rigid body modification 
    # pdb = '/home/schowell/data/myData/manualStructures/gH5_xTwist/pdb/gH5x4_m14_n2n4_r24.pdb'
    ## angle between sequenctial NCP z-axes:  [ 107.99115398  107.37279955  107.69404907]
    ## angle between stacked NCP z-axes:      [ 12.62481159  17.77834797]
    ## opening angles between NCPs:           [ 22.91152629  23.01043624]
    ## stack distances:                       [57.809280925176395, 58.203100113549787]
    
    # the twist starting structure
    # pdb = '/home/schowell/data/myData/manualStructures/gH5_xTwist/pdb/gH5x4_p0_n2n4_r24.pdb'
    ## angle between sequenctial NCP z-axes:  [ 94.00369639  93.60561237  93.97278084]
    ## angle between stacked NCP z-axes:      [ 12.62481159  13.19528317]
    ## opening angles between NCPs:           [ 22.91110313  23.01163443]
    ## stack distances:                       [57.809280925176395, 58.205898862508221]

    # get tetramer axes
    pkl_file = pdb[:-3] + 'pkl'
    try:
        pkl_in = open(pkl_file, 'rb')
        load_masks = True
        all_ncp_bases      = pickle.load(pkl_in)
        all_ncp_masks      = pickle.load(pkl_in)
        all_dyad_bases     = pickle.load(pkl_in)
        all_dyad_masks     = pickle.load(pkl_in)
        all_ref_masks      = pickle.load(pkl_in)
        all_ncp_origins    = pickle.load(pkl_in)
        all_ncp_axes       = pickle.load(pkl_in)
        all_ncp_opt_params = pickle.load(pkl_in)
        all_ncp_dyad_mol   = pickle.load(pkl_in)
        all_ncp_plot_vars  = pickle.load(pkl_in)
        pkl_in.close()
    except:
        all_ncp_plot_vars, all_ncp_axes, all_ncp_origins = get_tetramer_axes(
            pdb, ncp_dna_resids, dna_ids, ncp_dyad_resids, ncp_ref_atom_resids)

    # geometry.show_ncps(all_ncp_plot_vars)
    
    # Get the angle between x-axes
    all_ncp_axes_array = np.array(all_ncp_axes)
    angle_btwn_z_axes = geometry.angle_btwn_v1_v2(all_ncp_axes_array[:-1,2,:],
                                                  all_ncp_axes_array[1:,2,:])
    print 'angle between sequenctial NCP z-axes: ', angle_btwn_z_axes[0]
    angle_btwn_z_axes2 = geometry.angle_btwn_v1_v2(all_ncp_axes_array[:-2,2,:],
                                                  all_ncp_axes_array[2:,2,:])
    print 'angle between stacked NCP z-axes: ', angle_btwn_z_axes2[0]
 
    # Get the opening angles
    all_ncp_origins_array = np.array(all_ncp_origins)
    v_ncp1_to_ncp2 = all_ncp_origins_array[1:] - all_ncp_origins_array[:-1]
    opening_angles = geometry.angle_btwn_v1_v2(-v_ncp1_to_ncp2[:-1], 
                                               v_ncp1_to_ncp2[1:])
    print 'opening angles between NCPs: ', opening_angles[0]
    
    # Get the NCP stack distance
    stack_distance = all_ncp_origins_array[:2] - all_ncp_origins_array[2:] 
    stack_distance = [np.sqrt(vec.dot(vec)) for vec in stack_distance]
    print 'stack distances: ', stack_distance

    # Get the array center
    center = np.array(all_ncp_origins_array).mean(axis=0)
    print 'tetranucleosome center: ', center

    
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

    phi, dxyz, plot_title = get_tetramer_angles(all_ncp_axes, all_ncp_origins)
    # plot all NCPs:
    geometry.show_ncp_geometry(all_ncp_plot_vars)
    geometry.show_ncps(all_ncp_plot_vars, title=plot_title)

    # plot NCP1-NCP2:
    geometry.show_ncps(all_ncp_plot_vars[:2], title=plot_title)

    # plot NCP1-NCP3:
    geometry.show_ncps(all_ncp_plot_vars[:3], title=plot_title)

    n1_axes_name = pdb[:-4] + '_n1_axes.txt'
    np.savetxt(n1_axes_name, all_ncp_axes[0], fmt='%1.6e')
    
    n2_axes_name = pdb[:-4] + '_n2_axes.txt'
    np.savetxt(n2_axes_name, all_ncp_axes[1], fmt='%1.6e')    

    n3_axes_name = pdb[:-4] + '_n3_axes.txt'
    np.savetxt(n3_axes_name, all_ncp_axes[1], fmt='%1.6e')    

    n4_axes_name = pdb[:-4] + '_n4_axes.txt'
    np.savetxt(n4_axes_name, all_ncp_axes[1], fmt='%1.6e')    
    
    origins_name = pdb[:-4] + '_orig.txt'
    np.savetxt(origins_name, all_ncp_origins, fmt='%1.6e')
    
    dxyz_name = pdb[:-4] + '_dxyz.txt'
    np.savetxt(dxyz_name, dxyz, fmt='%1.6e')
    
    # r1_name = pdb[:-4] + '_r1.txt'
    # np.savetxt(r1_name, r1, fmt='%1.6e')

    # r2_name = pdb[:-4] + '_r2.txt'
    # np.savetxt(r2_name, r2, fmt='%1.6e')
    
    # h_name = pdb[:-4] + '_h.txt'
    # np.savetxt(h_name, h, fmt='%1.6e')

    phi_name = pdb[:-4] + '_phi.txt'
    np.savetxt(phi_name, phi, fmt='%0.3f')
    
    # data = 'N4merH5TE_zeroCon.iq'
    # data = 'N4merH5TE_zeroCon.i0q'
    # data = 'N4merH5Mg1_zeroCon.iq'
    # data = 'N4merH5Mg1_zeroCon.i0q'
    print '\m/ >.< \m/'    
    