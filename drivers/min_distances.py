#!/usr/bin/env python
#!/share/apps/bin/python
#
# Author:  Steven C. Howell
# Purpose: Generate a histogram of distances from pdb files
# Created: 12 May 2015
#
# $Id: dna_mc_driver.py 128 2015-02-23 18:52:00Z xraylab $
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

# def not_hydrogen(mol):
    # i = 0
    # for i, element in enumerate(mol.elements()):
        # if i == 0:
            # basis = ''        
        # else:
            # basis += ' or '
        # if element is not 'H':
            # basis += 'element %s' % element
    # return basis

import sassie.sasmol.sasmol as sasmol
import x_dna.energy.collision as collision
import numpy as np
import x_dna.util.basis_to_python as b2p
import matplotlib.pyplot as plt
import os.path as op
# import h5py
import pandas as pd

# f_name = ['new_dsDNA60.pdb'] 
f_name = ['/home/schowell/data/myData/sassieRuns/2x167/dimer_mod.pdb'] 
# f_name = ['/home/schowell/data/myData/sassieRuns/2x167/dimer_mod_d.pdb'] 
# f_name = ['/home/schowell/data/myData/sassieRuns/3x167/3x167.pdb'] 
# f_name = ['/home/schowell/data/myData/sassieRuns/4x167/distances/c11a.pdb']
# f_name = ['/home/schowell/data/myData/sassieRuns/4x167/distances/c11b.pdb']
# f_name = ['/home/schowell/data/myData/sassieRuns/4x167/c11.pdb']
# f_name = ['/home/schowell/data/myData/1kx5/1kx5_ncp.pdb'] 
# f_name = ['/home/schowell/data/myData/1zbb/1ZBB.pdb']

# minimized pdb files:
# f_name = ['/home/schowell/data/code/pylib/x_dna/build_mol/c11_v2/c11_val/'
          # 'minimization/c11_min.pdb']

# N_closest = 100
cutoff = 1.2

for pdb in f_name:
    in_mol = sasmol.SasMol(0)
    in_mol.read_pdb(pdb)

    #s n_atoms = mol.natoms()
    #s coor = mol.coor()
    #s dist = np.zeros((n_atoms, n_atoms))
    #s all_dist_m = collision.distances(coor, dist)
    #s all_dist = all_dist_m[all_dist_m>0]
    #s all_min = all_dist.min()
    #s n0, bins0, patches0 = plt.hist(all_dist, 500, alpha=0.75, facecolor='green',
    #s                                label='min all: %f' % all_min)    
   
    freq_file = '%s_freq%s.txt' % (pdb[:-4], str(cutoff))
    basis_filter = "name[i][0] != 'H' "
    e, mask = in_mol.get_subset_mask(basis_filter)
    mol = sasmol.SasMol(0)
    e = in_mol.copy_molecule_using_mask(mol, mask, 0)
    if op.exists(freq_file) and False:
        freq_df = pd.read_csv(freq_file, sep='\t')
    else:
        coor = np.array(mol.coor(), dtype=np.float32, order='F')
        n_heavy = coor.shape[1]
        # n_dist = n_heavy**2 #
        n_dist = n_heavy*(n_heavy-1)/2
        heavy_dist = np.zeros((n_heavy, n_heavy), dtype=np.float32, order='F')
        n_coor = n_heavy*3
        n_floats = n_dist + n_coor
        n_GB = n_floats * 4 / float(2**30) *2 # *2 b/c it will use this much mem
        n_Gb = n_floats * 4 / 1E9 *2          # in both python and fortran
        print 'calculating %d distances between %d atoms' % (n_dist, n_heavy)
        print 'should use between %f and %f GB of memory' % (n_GB, n_Gb)
        heavy_dist = collision.distances(coor, heavy_dist)
        # h5f = h5py.File(dist_file, 'w')
        # h5f.create_dataset('distances', data=heavy_dist)
        # h5f.close()
    
        # too much memory
        # max_dist = heavy_dist.max()
        # heavy_dist[np.tril_indices(n_heavy)] = max_dist * 2 
        # i_sort = heavy_dist.argsort(axis=None)[:N_closest]
        # r, c = np.unravel_index(i_sort, heavy_dist.shape)
        # del i_sort
        
        # much less memory
        print 'creating mask'
        mask = (heavy_dist > 0) * (heavy_dist < cutoff)
        
        print 'finding coordinates'
        r, c = np.where(mask)
        del mask
        
        dist = heavy_dist[r, c]
        print 'removing zeros'
        heavy_dist = heavy_dist[heavy_dist > 0] # linearizes the array
        # heavy_dist = heavy_dist[heavy_dist<=max_dist]
        
        print 'creating DataFrames'
        # i_dist = np.triu_indices(n_heavy,1)
        seg1  = [mol.segname()[i] for i in r]
        name1 = [mol.name()[i]    for i in r]
        res1  = [mol.resid()[i]   for i in r]
        seg2  = [mol.segname()[i] for i in c]
        name2 = [mol.name()[i]    for i in c]
        res2  = [mol.resid()[i]   for i in c]
        ni1 = [sum(r[i] == r) for i in xrange(len(r))]
        ni2 = [sum(c[i] == c) for i in xrange(len(c))]
        
        dist_dict = {'i1': r, 'i2': c, 'distance': dist, 'segname1': seg1, 
                     'segname2': seg2, 'resid1': res1, 'resid2': res2, 'name1': 
                     name1, 'name2': name2, 'n i1': ni1, 'n i2': ni2}
        dist_df = pd.DataFrame(dist_dict)
        dist_df.sort('distance', inplace=True)
        dist_file = '%s_dist%s.txt' % (pdb[:-4], str(cutoff))
        dist_df.to_csv(dist_file, sep='\t')
    
        # Get the most frequent close atoms
        i_close = np.concatenate((r, c))
        i_count = np.bincount(i_close)
        ii = np.nonzero(i_count)[0]
        i_binned = np.vstack((ii,i_count[ii])).T
        i_binned[i_binned[:,1].argsort()[::-1]]
        i_atom = i_binned[:,0]
        atom_frequency = i_binned[:,1]
        
        segname = [mol.segname()[i] for i in i_atom]
        name    = [mol.name()[i] for i in i_atom]
        resid   = [mol.resid()[i] for i in i_atom]
        mn_dist = [np.nanmin([dist_df.loc[dist_df['i2'] == i]['distance'].max(),
                              dist_df.loc[dist_df['i1'] == i]['distance'].max()]
                             ) for i in i_atom]
        in_atom1 = [sum(dist_df['i1'] == i_atom[i]) for i in 
                    xrange(len(i_atom))]
        in_atom2 = [sum(dist_df['i2'] == i_atom[i]) for i in 
                    xrange(len(i_atom))]

        freq_dict = {'i': i_atom, 'frequency': atom_frequency, 'resid': resid,
                      'segname': segname, 'name': name, 'distance': mn_dist,
                      'in atom1': in_atom1, 'in atom2': in_atom2}
        freq_df = pd.DataFrame(freq_dict)
        freq_df.sort(['frequency', 'distance'], inplace=True, 
                     ascending=[False, True])
        freq_df.to_csv(freq_file, sep='\t')
    
    
        print 'getting min'
        heavy_min = heavy_dist.min()
    
        print 'subplot 1'
        plt.subplot(211)
        n1, bins1, patches1 = plt.hist(heavy_dist, 50, alpha=0.75, 
                                       facecolor='blue', label='heavy atoms' %
                                       heavy_min)
        plt.ylabel('Counts')
        plt.title(r'Atomic distances: %s' % op.split(pdb)[1])
        plt.grid(True)
        lg = plt.legend()
        lg.draw_frame(False)
        
        print 'subplot 2'
        plt.subplot(212)
        # n1, bins1, patches1 = plt.hist(heavy_dist[heavy_dist<3], 
        n1, bins1, patches1 = plt.hist(heavy_dist[np.argpartition(
            heavy_dist, 10)[:10]], 5, alpha=0.75, facecolor='green', 
                                       label='min heavy: %f' % heavy_min)
        plt.xlabel('Distances')
        plt.ylabel('Counts')
        plt.grid(True)
        lg = plt.legend(loc=9)
        lg.draw_frame(False)
    
        print 'saving figures'
        plt.savefig(pdb[:-3] + 'png', dpi=400, bbox_inches='tight')
        plt.savefig(pdb[:-3] + 'eps', bbox_inches='tight')    
        # plt.show()
        
    print 'creating pdb without these atoms'
    # this will leave out all the hydrogen atoms (faster data transfer)
    e, mask = mol.get_subset_mask('all')
    for i in xrange(len(dist_df)):
        if dist_df['name1'][i] in ['CA', 'P']:
            mask[dist_df['i2'][i]] = 0
        elif dist_df['name2'][i] in ['CA', 'P']:
            mask[dist_df['i1'][i]] = 0
        else:
            if dist_df['n i2'][i] > dist_df['n i1'][i]:
                mask[dist_df['i2'][i]] = 0
            else:
                mask[dist_df['i1'][i]] = 0

    new_mol = sasmol.SasMol(0)
    error = mol.copy_molecule_using_mask(new_mol, mask, 0)
    out_pdb = '%s_d%s.pdb' % (pdb[:-4], str(cutoff))
    new_mol.write_pdb(out_pdb, 0, 'w')
    new_mol.read_pdb(out_pdb)   
    new_mol.write_pdb(out_pdb, 0, 'w')
    
print '\m/ >.< \m/'
