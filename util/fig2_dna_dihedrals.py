#!/usr/bin/python
# Author:  Steven C. Howell
# Purpose: script for calculating DNA dihedral angles and generating plots
# Created: 22 September 2014
# $Id: $

import sys
import os
import logging
import subprocess
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import sassie.calculate.dna_dihedral as dd
import x_dna.util.gw_plot as gwp
import multiprocessing as mp

def main():
    if '-v' in sys.argv:
        logging.basicConfig(filename='%s.log' % __file__[:-3],
                            level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
        # print 'logging level set to DEBUG'
    else:
        logging.basicConfig()

    good_dihedral = False
    scatter_dihedrals = False
    calc_dihedrals = False
    parallel = False
    frequency = -1

    # ~~~~~~~~ FILE INPUT ~~~~~~~~~~ #
    # first and last DNA resdue of one of the DNA chains
    # first_last_resids = [[15, 161], [339, 193]]
    # run_dir = '/home/schowell/data/myData/my_mono/mono_ncp_gb1/energy_minimization/'
    # pdb_file_name = run_dir + 'min_run_0.dcd.pdb'
    # flex_file = run_dir + 'min_run_0_arm1.flex'
    # flex_file = run_dir + 'min_run_0_arm2.flex'
    # drude = False
    # dcd_files = []
    # dcd_files.append(run_dir + 'min_run_0.dcd')


    # first_last_resids = [[15, 161], [339, 193]]
    # run_dir = '/home/schowell/data/myData/my_mono/mono_ncp_gb50/energy_minimization/'
    # pdb_file_name = run_dir + 'min_run_0.dcd.pdb'
    # flex_file = run_dir + 'min_run_0_arm1.flex'
    # flex_file = run_dir + 'min_run_0_arm2.flex'
    # drude = False
    # dcd_files = []
    # dcd_files.append(run_dir + 'min_run_0.dcd')

    #~~~ new charmm36 60 bps dsDNA (compare before and after min) ~~~#
    # first_last_resids = [[1, 60], [61, 120]]
    # run_dir = '/home/schowell/data/code/pylib/sassie_2_na/build_mol/dsDNA60/'
    # flex_file = '/home/schowell/data/myData/dihedrals/dsDNA_60bps/new_dsDNA60.flex'
    # drude = False
    # dcd_files = []
    # pdb_file_name = run_dir + 'c36_dna_raw.pdb'
    # dcd_files.append(run_dir + 'c36_dna_raw_x2.dcd')
    # pdb_file_name = run_dir + 'output_building/c36_dna60.pdb'
    # dcd_files.append(run_dir + 'output_building/c36_dna60_raw_min.dcd')

    #~~~ select scatter plot ~~~#
    first_last_resids = [[1, 60], [61, 120]]
    run_dir = '/home/schowell/data/myData/dihedrals/dsDNA_60bps/'
    pdb_file_name = run_dir + 'c36_min_dsDNA60.pdb'
    flex_file = run_dir + 'new_dsDNA60.flex'
    drude = False
    dcd_files = [run_dir + 'ma10_mc1_1000spb_715_raw_min.dcd']
    prefix = 'raw_min_MD'

    #~~~ determine how many time DNA can be minimized ~~~#
    # first_last_resids = [[1, 60], [61, 120]]
    # run_dir = '/home/schowell/data/myData/dihedrals/dsDNA_60bps/'
    # pdb_file_name = run_dir + 'c36_min_dsDNA60.pdb'
    # flex_file = run_dir + 'new_dsDNA60.flex'
    # drude = False
    # dcd_files = []
    # dcd_files.append(run_dir + 'ma0_spb1/no_mc/mmm_x100.dcd')
    # dcd_files.append(run_dir + 'm_ma0_spb1/no_mc/m_x100.dcd')

    #~~~ 60 bps dsDNA ~~~#
    # first_last_resids = [[1, 60], [61, 120]]
    # run_dir = '/home/schowell/data/myData/dihedrals/dsDNA_60bps/'
    # pdb_file_name = run_dir + 'c36_min_dsDNA60.pdb'
    # flex_file = run_dir + 'new_dsDNA60.flex'
    # drude = False
    # dcd_files = []

    # dcd_files.append(run_dir + 'ma10_spb1000/mc1_spb1000_sparse.dcd')
    # dcd_files.append(run_dir + 'ma1_spb1000/mc1_spb1000_sparse.dcd')
    # dcd_files.append(run_dir + 'ma20_spb1000/mc1_spb1000_sparse.dcd')
    # dcd_files.append(run_dir + 'ma30_spb1000/mc1_spb1000_sparse.dcd')
    # dcd_files.append(run_dir + 'ma40_spb1000/mc1_spb1000_sparse.dcd')
    # dcd_files.append(run_dir + 'ma50_spb1000/mc1_spb1000_sparse.dcd')
    # dcd_files.append(run_dir + 'ma5_spb1000/mc1_spb1000_sparse.dcd')
    # dcd_files.append(run_dir + 'ma60_spb1000/mc1_spb1000_sparse.dcd')

    # first_last_resids = [[1, 60], [61, 120]]
    # run_dir = '/home/schowell/data/myData/dihedrals/dsDNA_60bps'
    # pdb_file_name = os.path.join(run_dir, 'c36_min_dsDNA60.pdb')
    # flex_file = os.path.join(run_dir, 'new_dsDNA60.flex')
    # drude = False
    # dcd_files = []

    # prefix = 'debug_MD'
    # run_dir = os.path.join(run_dir, 'ma20_spb1000')
    # dcd_files.append(os.path.join(run_dir, 'mc1_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'mc2_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'mmm_energy_minimization/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'mmm_energy_minimization_BK/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'm_energy_minimization/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'm_energy_minimization_BK/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'm_energy_minimization_redo/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'mmm_energy_minimization_redo/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'm_energy_minimization_v2/mc2_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'mmm_energy_minimization_v2/mc1_spb1000_sparse_mmm_tmp.dcd'))

    # prefix = 'raw_MD'
    # dcd_files.append(os.path.join(run_dir, 'ma60_spb1000/mc1_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma50_spb1000/mc1_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma40_spb1000/mc1_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma30_spb1000/mc1_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma20_spb1000/mc2_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma10_spb1000/mc1_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma5_spb1000/mc1_spb1000_sparse.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma1_spb1000/mc1_spb1000_sparse.dcd'))

    # prefix = 'mmm_MD'
    # dcd_files.append(os.path.join(run_dir, 'ma60_spb1000/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma50_spb1000/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma40_spb1000/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma30_spb1000/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma20_spb1000/mc2_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma10_spb1000/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma5_spb1000/mc1_spb1000_sparse_mmm.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma1_spb1000/mc1_spb1000_sparse_mmm.dcd'))

    # prefix = 'm_MD'
    # dcd_files.append(os.path.join(run_dir, 'ma60_spb1000/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma50_spb1000/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma40_spb1000/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma30_spb1000/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma20_spb1000/mc2_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma10_spb1000/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma5_spb1000/mc1_spb1000_sparse_m.dcd'))
    # dcd_files.append(os.path.join(run_dir, 'ma1_spb1000/mc1_spb1000_sparse_m.dcd'))

    # prefix = 'sanity_MD'
    # dcd_files.append(os.path.join(run_dir, 'ma20_spb1000/sanity/monte_carlo/'
                                  # 'sanity_ma20_1k.dcd'))

    # run_dir = 'charmm36/'
    # first_last_resids = [[1,12],[13,24]]

    # drude = True
    # dcd_files = []
    # flex_file = run_dir + 'test_dna.flex'
    # pdb_file_name = run_dir + 'test_dna.pdb'
    # dcd_files.append(run_dir + 'drude/step5_equilibration.dcd')
    # # dcd_files.append(run_dir + 'step5_equil_short.dcd')
    # # dcd_files.append(run_dir + 'step5_equilibration.dcd')

    # drude = False
    # dcd_files = []
    # flex_file = run_dir + '/non_polarized/dyn1.flex'
    # pdb_file_name = run_dir + '/non_polarized/dyn1.pdb'
    # dcd_files.append(run_dir + 'dyn2.dcd')

    # run_dir = 'test_dna/'
    # first_last_resids = [[1, 12], [13, 24]]
    # pdb_file_name = run_dir + 'new_dna.pdb'
    # flex_file = run_dir + 'dna.flex'
    # dcd_files.append(run_dir + 'raw_min.dcd')

    # first_last_resids = [[1, 60], [61, 120]]
    # run_dir = 'test_dna/raw_version/'
    # pdb_file_name =  run_dir + 'new_dna.pdb'
    # flex_file = run_dir + 'new_dna.flex'
    # dcd_files = []
    # dcd_files.append(run_dir + 'output/dyn1_dna.dcd')
    # ~~~~~~~~ END FILE INPUT ~~~~~~~~~~ #

    # ~~~~~~~~ RUN INPUT ~~~~~~~~~~ #
    calc_dihedrals = True
    good_dihedral = True
    scatter_dihedrals = True
    # parallel = True
    # ~~~~~~~~ RUN INPUT ~~~~~~~~~~ #

    processes = []
    used_goback = False
    max_frame = None
    scale_spb = 100.0 / 58.0
    # scale_spb = 1.0
    show = False

    for dcd_file_name in dcd_files:
        # calculate and save the DNA dihedral angels
        if calc_dihedrals:
            hdf_file = dcd_file_name[:-3] + 'hdf'
            if os.path.exists(hdf_file):
                print "HDF file already exists: %s" % hdf_file
                print "skipping calculations (remove file to recalculate)"
            else:
                print "calculating dihedral angles for %s" % hdf_file
                if parallel:
                    get_angles_inputs = (dcd_file_name, pdb_file_name,
                                         first_last_resids, flex_file)

                    processes.append(mp.Process(target=dd.get_angles_df,
                                                args=get_angles_inputs))
                else:
                    angles_df = dd.get_angles_df(dcd_file_name, pdb_file_name,
                                                 first_last_resids, flex_file,
                                                 drude=drude)
                    print "finished angle calculations for %s" % hdf_file

        if good_dihedral:
            pub_scatter_dihedrals(dcd_file_name, frequency, scale_spb, True)
            if parallel:
                plot_inputs = (dcd_file_name, max_frame, scale_spb,
                               used_goback, show, True)
                processes.append(mp.Process(target=plot_good_dihedrals,
                                            args=plot_inputs))
                # plot_inputs = (dcd_file_name, max_frame, scale_spb,
                               # used_goback, show, False)
                # processes.append(mp.Process(target=plot_good_dihedrals,
                                            # args=plot_inputs))
            else:
                plot_good_dihedrals(dcd_file_name, max_frame, scale_spb,
                                    used_goback, show, True)
                # plot_good_dihedrals(dcd_file_name, max_frame, scale_spb,
                                    # used_goback, show, False)

        if scatter_dihedrals:
            if parallel:
                scatter_inputs = (dcd_file_name, frequency, scale_spb, True)
                processes.append(mp.Process(target=scatter_plot_dihedrals,
                                          args=scatter_inputs))
                scatter_inputs = (dcd_file_name, frequency, scale_spb, False)
                processes.append(mp.Process(target=scatter_plot_dihedrals,
                                          args=scatter_inputs))
            else:
                scatter_plot_dihedrals(dcd_file_name, frequency, scale_spb,
                                       True)
                scatter_plot_dihedrals(dcd_file_name, frequency, scale_spb,
                                       False)


    if parallel:
        for p in processes:
            p.start()

        for p in processes:
            p.join()

    # compare_average_in_range(dcd_files)

def plot_good_dihedrals(dcd_file_name, max_frame=None, scale_spb=1,
                        used_goback=False, show=False, MD=True):
    if MD:
        if max_frame:
            npy_file = dcd_file_name[:-4] + '_MD_%d.npy' % max_frame
        else:
            npy_file = dcd_file_name[:-4] + '_MD.npy'
    else:
        if max_frame:
            npy_file = dcd_file_name[:-4] + '_%d.npy' % max_frame
        else:
            npy_file = dcd_file_name[:-4] + '.npy'

    try:
        angle_vs_steps = np.load(npy_file)
        print 'loaded percent of good angles file: %s' % npy_file

    except:
        print 'failed to load percent of good angles file: %s' % npy_file

        angles = pd.read_hdf(dcd_file_name[:-3] + 'hdf', 'angles')
        angles['-90'] = np.abs(-90 - (angles['epsilon'] - angles['zeta']))
        angles['+90'] = np.abs(90 - (angles['epsilon'] - angles['zeta']))
        angles['min'] = angles.loc[:,['+90', '-90']].min(axis=1)
        angles['B-type'] = 'B'
        angles['B-type'][angles['min'] == angles['+90']] = 'BII'
        angles['B-type'][angles['min'] == angles['-90']] = 'BI'

        try:
            scale_spb == 1  # enforce there is not a scale
            all_spb = np.loadtxt(dcd_file_name[:-3] + 'spb')
            print ('loaded steps/bead file: %s' %
                   (dcd_file_name[:-3] + 'spb'))

        except:
            print ('failed to load steps/bead file: %s \n'
                   'using %f step/s between frames' %
                   (dcd_file_name[:-3] + 'spb', scale_spb))
            all_spb = np.unique(angles['frame'])
            all_spb = all_spb.reshape(len(all_spb), 1)
            all_spb = np.concatenate((all_spb, all_spb*scale_spb), axis=1)

        if max_frame:
            nf = max_frame
        else:
            nf = len(all_spb)
        spb, indices = np.unique(all_spb[:nf, 1], return_inverse=True)
        n_steps = len(spb)

        angle_vs_steps = np.zeros((n_steps, 8))  # 7 angles
        angle_vs_steps[:, 0] = spb

        if MD:
            limit_file = ('/home/schowell/data/myData/dihedrals/'
                          'dsDNA_60bps/md_angles.txt')
        else:
            limit_file = ('/home/schowell/data/myData/dihedrals/'
                          'dsDNA_60bps/dna_angles.txt')
        limits = pd.read_csv(limit_file, index_col=0)
        limits['low'] = limits['mean'] - 2 * limits['sd']
        limits['high'] = limits['mean'] + 2 * limits['sd']

        for (i, steps) in enumerate(spb):
            if used_goback:
                df, _ = get_step_i_dihedrals(i, indices, angles)
            else:
                df = angles[angles['frame'] == i+1]

            beta_good = sum(
                ((df['beta'] > limits.loc['beta1']['low']) &
                 (df['beta'] < limits.loc['beta1']['high']) &
                 (df['B-type'] == 'BI')) |
                ((df['beta'] > limits.loc['beta2']['low']) &
                 (df['beta'] < limits.loc['beta2']['high']) &
                 (df['B-type'] == 'BII')))
            n_beta = df['beta'].notnull().sum() * 1.0
            angle_vs_steps[i, 2] = beta_good / n_beta

            gamma_good = sum(
                (df['gamma'] > limits.loc['gamma']['low']) &
                (df['gamma'] < limits.loc['gamma']['high']))
            n_gamma = df['gamma'].notnull().sum() * 1.0
            angle_vs_steps[i, 3] = gamma_good / n_gamma

            delta_good = sum(
                ((df['delta'] > limits.loc['delta1']['low']) &
                 (df['delta'] < limits.loc['delta1']['high']) &
                (df['B-type'] == 'BI'))  |
                ((df['delta'] > limits.loc['delta2']['low']) &
                 (df['delta'] < limits.loc['delta2']['high']) &
                 (df['B-type'] == 'BII')))
            n_delta = df['delta'].notnull().sum() * 1.0
            angle_vs_steps[i, 4] = delta_good / n_delta

            chi_good = sum(
                (((df['resname'] == 'C') | (df['resname'] == 'T')) &
                 (df['chi'] > limits.loc['chi1py']['low']) &
                 (df['chi'] < limits.loc['chi1py']['high']) &
                 (df['B-type'] == 'BI'))  |
                (((df['resname'] == 'G') | (df['resname'] == 'A')) &
                 (df['chi'] > limits.loc['chi1pu']['low']) &
                 (df['chi'] < limits.loc['chi1pu']['high']) &
                 (df['B-type'] == 'BI'))  |
                ((df['chi'] > limits.loc['chi2']['low']) &
                 (df['chi'] < limits.loc['chi2']['high']) &
                 (df['B-type'] == 'BII')))
            n_chi = df['chi'].notnull().sum() * 1.0
            angle_vs_steps[i, 7] = chi_good / n_chi

            alpha_good = sum(
                (df['alpha'] > limits.loc['alpha']['low']) &
                (df['alpha'] < limits.loc['alpha']['high']))
            n_alpha = df['alpha'].notnull().sum() * 1.0
            angle_vs_steps[i, 1] = alpha_good / n_alpha

            epsilon_good = sum(
                ((df['epsilon'] > limits.loc['epsilon1']['low']) &
                 (df['epsilon'] < limits.loc['epsilon1']['high']) &
                 (df['B-type'] == 'BI'))  |
                ((df['epsilon'] > limits.loc['epsilon2']['low']) &
                 (df['epsilon'] < limits.loc['epsilon2']['high']) &
                 (df['B-type'] == 'BII')))
            n_epsilon = df['epsilon'].notnull().sum() * 1.0
            angle_vs_steps[i, 5] = epsilon_good / n_epsilon

            zeta_good = sum(
                ((df['zeta'] > limits.loc['zeta1']['low']) &
                 (df['zeta'] < limits.loc['zeta1']['high']) &
                 (df['B-type'] == 'BI'))  |
                ((df['zeta'] > limits.loc['zeta2']['low']) &
                 (df['zeta'] < limits.loc['zeta2']['high']) &
                 (df['B-type'] == 'BII')))
            n_zeta = df['zeta'].notnull().sum() * 1.0
            angle_vs_steps[i, 6] = zeta_good / n_zeta

        np.save(npy_file, angle_vs_steps)

    fig = plt.figure(figsize=(6, 4))
    gs1 = gridspec.GridSpec(1, 1, left=0.075, right=0.75, wspace=0.1,
                            hspace=0, top=0.95)
    ax1 = plt.subplot(gs1[:, 0])
    angle_labels = [r'$\alpha$', r'$\beta$', r'$\gamma$', r'$\delta$',
                    r'$\epsilon$', r'$\zeta$', r'$\chi$']
    plot_order = [2, 3, 4, 7, 1, 5, 6]
    # plot_order = [1, 5, 6]
    # plot_order = [2, 3, 4, 7]
    ax1.plot([0, 10000], [95, 95], '--', c='DimGray', linewidth=2)
             # label=r'$2\sigma$ cutoff')
    if len(angle_vs_steps)>2:
        for (c, i_angle) in enumerate(plot_order):
            ax1.plot(angle_vs_steps[:, 0], angle_vs_steps[:, i_angle]*100, '-',
                     c=gwp.qual_color(c), label=angle_labels[i_angle - 1])
        ax1.set_ylim([73, 102])
        ax1.set_xscale('log')
        x_max = angle_vs_steps[-1, 0]
        ax1.set_xlim([0, x_max])
    else:
        for (c, i_angle) in enumerate(plot_order):
            ax1.plot(angle_vs_steps[:, 0], angle_vs_steps[:, i_angle]*100, '-',
                     color=gwp.qual_color(c), marker=gwp.symbol_order(c),
                     mec=gwp.qual_color(c), label=angle_labels[i_angle - 1],
                     markersize=15, mfc='none')
        ax1.set_ylim([-1, 110])

    ax1.set_ylabel(r'% in range')
    ax1.set_xlabel(r'Minimization Iterations')
    default_fontsize = 12
    lg = plt.legend(loc='lower left', scatterpoints=1, numpoints=1,
                    prop={'size': default_fontsize})
    lg.draw_frame(False)
    if show:
        plt.show()
    else:
        save_name = '%s_in_range' % dcd_file_name[:-4]
        if MD:
            save_name += '_MD'
        if max_frame:
            save_name += '_' + str(max_frame)
        print 'View figure: \nevince %s &' % save_name
        plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
        plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')

    print 'done plotting'

    all_in_range_file = dcd_file_name[:-4] + '_ave_in_range'
    if MD:
        all_in_range_file += '_MD'
    all_in_range = np.copy(angle_vs_steps[:,:2])
    all_in_range[:,1] = angle_vs_steps[:,1:].mean(axis=1)
    np.savetxt(all_in_range_file + '.out', all_in_range)


def compare_average_in_range(dcd_files, MD=True, show=False):

    fig = plt.figure(figsize=(6, 4))
    gs1 = gridspec.GridSpec(1, 1, left=0.075, right=0.75, wspace=0.1,
                            hspace=0, top=0.95)
    ax1 = plt.subplot(gs1[:, 0])
    labels = [r'$\delta\theta=60$', r'$\delta\theta=50$', r'$\delta\theta=40$',
              r'$\delta\theta=30$', r'$\delta\theta=20$', r'$\delta\theta=10$',
              r'$\delta\theta=5$', r'$\delta\theta=1$']

    plt.axhline(95, ls='--', c='DimGray', linewidth=2)
    all_in_range = []
    max_steps = []
    min_steps = []
    for (i, dcd_file) in enumerate(dcd_files):
        all_in_range_file = os.path.join(os.path.split(dcd_file)[0],
                                    'average_in_range')
        if MD:
            all_in_range_file += '_MD'
        all_in_range.append(np.loadtxt(all_in_range_file + '.out'))
        max_steps.append(all_in_range[i][-1, 0])
        min_steps.append(all_in_range[i][0, 0])

        ax1.plot(all_in_range[i][:, 0], all_in_range[i][:, 1]*100, '-',
                 c=gwp.qual_color(i), label=labels[i])

    ax1.set_xscale('log')
    x_min = min(min_steps)
    x_max = min(max_steps)
    ax1.set_xlim([x_min, x_max])
    ax1.set_ylim([73, 102])

    ax1.set_ylabel(r'% in Range')
    ax1.set_xlabel(r'Minimization Iterations')
    default_fontsize = 12
    lg = plt.legend(loc='lower left', scatterpoints=1, numpoints=1,
                    prop={'size': default_fontsize})
    lg.draw_frame(False)
    if show:
        plt.show()
    else:
        save_name = '/home/schowell/data/myData/dihedrals/dsDNA_60bps/average_in_range'
        if MD:
            save_name += '_MD'
        print 'View figure: \nevince %s &' % save_name
        plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
        plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')

def get_step_i_dihedrals(i, indices, all_df):
    df_list = []
    frames = np.where(indices == i)[0]
    for i_angle in frames:
        df_list.append(all_df[all_df['frame'] == i_angle + 1])
    i_df = pd.concat(df_list)
    return i_df, len(frames)


def limit_patch(x_key, y_key, limits, ax, b1=1):
    red_alpha_p5 = np.array([235, 157, 158], dtype=float) / 255
    yellow_alpha_p5 = np.array([241, 218, 151], dtype=float) / 255
    overlap_red_yellow = np.array([228, 182, 48], dtype=float) / 255
    if b1 == 1:
        color = red_alpha_p5
    elif b1 == 2 :
        color = yellow_alpha_p5
    elif b1 == 3:
        color = overlap_red_yellow
    limit_patch = patches.Rectangle(
        (limits.loc[x_key]['low'], limits.loc[y_key]['low']),
        4 * limits.loc[x_key]['sd'], 4 * limits.loc[y_key]['sd'],
        edgecolor='none', alpha=1, facecolor=color)
    ax.add_patch(limit_patch)



def pub_scatter_dihedrals(dcd_file_name, frequency=-1, scale_spb=1, MD=True):
    prefix = dcd_file_name[:-4] + '_'
    angles = pd.read_hdf(dcd_file_name[:-3] + 'hdf', 'angles')

    try:
        all_spb = np.loadtxt(dcd_file_name[:-3] + 'spb')
    except:
        all_spb = np.unique(angles['frame'])
        all_spb = all_spb.reshape(len(all_spb), 1)
        all_spb = np.concatenate((all_spb, all_spb*scale_spb), axis=1)

    spb, indices = np.unique(all_spb[:, 1], return_inverse=True)
    n_steps = len(spb)

    if MD:
        limit_file = ('/home/schowell/data/myData/dihedrals/'
                      'dsDNA_60bps/md_angles.txt')
    else:
        limit_file = ('/home/schowell/data/myData/dihedrals/'
                      'dsDNA_60bps/dna_angles.txt')
    limits = pd.read_csv(limit_file, index_col=0)
    limits['low'] = limits['mean'] - 2 * limits['sd']
    limits['high'] = limits['mean'] + 2 * limits['sd']

    angle_label = {'alpha': r'$\alpha$', 'beta': r'$\beta$',
                   'gamma': r'$\gamma$', 'delta': r'$\delta$',
                   'epsilon': r'$\epsilon$', 'zeta': r'$\zeta$',
                   'chi': r'$\chi$'}
    # matplotlib.rcParams['xtick.direction'] = 'in'
    # matplotlib.rcParams['ytick.direction'] = 'in'

    if n_steps == 2:
        df_raw, _ = get_step_i_dihedrals(0, indices, angles)
        ax_array = make_pub_scatter_plots(limits, df_raw, angle_label)

        df_min, _ = get_step_i_dihedrals(1, indices, angles)
        ax_array = make_pub_scatter_plots(limits, df_min, angle_label,
                                          i_color=2, symbol='s',
                                          ax_array=ax_array)

        ax_array[1].legend(['Raw', 'Minimized'], loc='lower right', numpoints=1)

        # plt.suptitle('Comparative scatter plot of selected torsional angles')
        save_name = prefix + 'dihedrals_pub'
        if MD:
            save_name += '_MD'
        plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
        plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
        plt.show()
    else:
        if frequency == -1:
            frequency = n_steps / 10
        for step in xrange(frequency-1, n_steps, frequency):
            df, nf = get_step_i_dihedrals(step, indices, angles)
            make_pub_scatter_plots(limits, df, angle_label)

            plt.suptitle('%d Steps (%d frame/s): scatter plots of selected'
                         'torsional angles' % (step + 1, nf))
            save_name = prefix + 'dihedrals_%dsteps' % (step + 1)
            if MD:
                save_name += '_MD'
            plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
            plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
        plt.show()
    print 'View publication scatter plot: \nevince %s.eps &' % save_name
    print 'pause'




def scatter_plot_dihedrals(dcd_file_name, frequency=-1, scale_spb=1, MD=True):
    prefix = dcd_file_name[:-4] + '_'
    angles = pd.read_hdf(dcd_file_name[:-3] + 'hdf', 'angles')

    try:
        all_spb = np.loadtxt(dcd_file_name[:-3] + 'spb')
    except:
        all_spb = np.unique(angles['frame'])
        all_spb = all_spb.reshape(len(all_spb), 1)
        all_spb = np.concatenate((all_spb, all_spb*scale_spb), axis=1)

    spb, indices = np.unique(all_spb[:, 1], return_inverse=True)
    n_steps = len(spb)

    if MD:
        limit_file = ('/home/schowell/data/myData/dihedrals/'
                      'dsDNA_60bps/md_angles.txt')
    else:
        limit_file = ('/home/schowell/data/myData/dihedrals/'
                      'dsDNA_60bps/dna_angles.txt')
    limits = pd.read_csv(limit_file, index_col=0)
    limits['low'] = limits['mean'] - 2 * limits['sd']
    limits['high'] = limits['mean'] + 2 * limits['sd']

    # plt.ion()  # enable interactive plotting (does not wait at plt.show())

    angle_label = {'alpha': r'$\alpha$', 'beta': r'$\beta$',
                   'gamma': r'$\gamma$', 'delta': r'$\delta$',
                   'epsilon': r'$\epsilon$', 'zeta': r'$\zeta$',
                   'chi': r'$\chi$'}
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'

    if n_steps == 2:
        df_0, _ = get_step_i_dihedrals(0, indices, angles)
        ax_array = make_selected_scatter_plots(limits, df_0, angle_label,
                                               symbol='+')

        df_1, _ = get_step_i_dihedrals(1, indices, angles)
        make_selected_scatter_plots(limits, df_1, angle_label, i_color=2,
                                    ax_array=ax_array)

        plt.legend(['raw', 'min'], loc='upper left',
                   bbox_to_anchor=[1.1475, 1.07], numpoints=1)

        plt.suptitle('Comparative scatter plot of selected torsional angles')
        save_name = prefix + 'dihedrals_comparison'
        if MD:
            save_name += '_MD'
        plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
        plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
        # plt.show()
    else:
        if frequency == -1:
            frequency = n_steps / 10
        for step in xrange(frequency-1, n_steps, frequency):
            df, nf = get_step_i_dihedrals(step, indices, angles)
            make_selected_scatter_plots(limits, df, angle_label)

            plt.suptitle('%d Steps (%d frame/s): scatter plots of selected'
                         'torsional angles' % (step + 1, nf))
            save_name = prefix + 'dihedrals_%dsteps' % (step + 1)
            if MD:
                save_name += '_MD'
            plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
            plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')

    # plt.show()
    print 'pause'

def plt_format(ax):

    ax.set_xlim([0, 360])
    ax.set_ylim([0, 360])
    # ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    # ax.set_yticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 90, 180, 270, 360])


def make_pub_scatter_plots(limits, df, angle_label, ax_array=[],
                                symbol='o', i_color=0):

    if len(ax_array) == 0:
        fig = plt.figure(figsize=(3, 5))
        gs = gridspec.GridSpec(2, 1, left=0.1, right=0.9, wspace=0, hspace=0)
        ax_array = [plt.subplot(gs[0]), plt.subplot(gs[1])]

    ax = ax_array[0]
    x = 'zeta'
    y = 'epsilon'
    limit_patch('zeta1', 'epsilon1', limits, ax, b1=1)
    limit_patch('zeta2', 'epsilon2', limits, ax, b1=2)
    ax.plot(df[x], df[y], symbol, mec=gwp.qual_color(i_color),
            mfc='none')
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    ax.tick_params(labelbottom='off')# ,  labeltop='on', labelright='on')

    ax = ax_array[1]
    x = 'zeta'
    y = 'alpha'
    limit_patch('zeta1', 'alpha', limits, ax, b1=1)
    limit_patch('zeta2', 'alpha', limits, ax, b1=2)
    ax.plot(df[x][0:-1], df[y][1:], symbol, mec=gwp.qual_color(i_color),
            mfc='none')
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)

    return ax_array


def make_selected_scatter_plots(limits, df, angle_label, ax_array=[],
                                symbol='x', i_color=0):

    if len(ax_array) == 0:
        fig, ax_array = plt.subplots(3, 3)

    ax = ax_array[0, 0]
    x = 'zeta'
    y = 'alpha'
    limit_patch('zeta1', 'alpha', limits, ax, b1=1)
    limit_patch('zeta2', 'alpha', limits, ax, b1=2)
    ax.plot(df[x][0:-1], df[y][1:], symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)
    ax.tick_params(labelbottom='off') # , labeltop='on')

    ax = ax_array[0, 1]
    x = 'zeta'
    y = 'beta'
    limit_patch('zeta1', 'beta1', limits, ax, b1=1)
    limit_patch('zeta1', 'beta2', limits, ax, b1=3)
    limit_patch('zeta2', 'beta1', limits, ax, b1=3)
    limit_patch('zeta2', 'beta2', limits, ax, b1=2)
    ax.plot(df[x][0:-1], df[y][1:], symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)
    ax.tick_params(labelleft='off', labelbottom='off') # , labeltop='on')

    ax = ax_array[0, 2]
    x = 'zeta'
    y = 'epsilon'
    limit_patch('zeta1', 'epsilon1', limits, ax, b1=1)
    limit_patch('zeta2', 'epsilon2', limits, ax, b1=2)
    ax.plot(df[x], df[y], symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    ax.tick_params(
        labelleft='off', labelbottom='off')# ,  labeltop='on', labelright='on')

    ax = ax_array[1, 0]
    x = 'alpha'
    y = 'gamma'
    limit_patch('alpha', 'gamma', limits, ax, b1=1)
    ax.plot(df[x], df[y], symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    ax.tick_params(labelbottom='off')

    ax = ax_array[1, 1]
    x = 'zeta'
    y = 'chi'
    limit_patch('zeta1', 'chi1pu', limits, ax, b1=1)
    limit_patch('zeta1', 'chi1py', limits, ax, b1=1)
    limit_patch('zeta2', 'chi2', limits, ax, b1=2)
    ax.plot(df[x], df[y],  symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    ax.tick_params(labelleft='off', labelbottom='off')

    ax = ax_array[1, 2]
    x = 'delta'
    y = 'chi'
    limit_patch('delta1', 'chi1pu', limits, ax, b1=1)
    limit_patch('delta1', 'chi1py', limits, ax, b1=1)
    limit_patch('delta2', 'chi2', limits, ax, b1=2)
    ax.plot(df[x], df[y],  symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    # , labelright='on')
    ax.tick_params(labelleft='off', labelbottom='off')

    ax = ax_array[2, 0]
    x = 'zeta'
    y = 'zeta'
    limit_patch('zeta1', 'zeta1', limits, ax, b1=1)
    limit_patch('zeta1', 'zeta2', limits, ax, b1=3)
    limit_patch('zeta2', 'zeta1', limits, ax, b1=3)
    # limit_patch('zeta2', 'zeta2', limits, ax, b1=True)
    ax.plot(df[x][:-1], df[y][1:],  symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)

    ax = ax_array[2, 1]
    x = 'epsilon'
    y = 'epsilon'
    limit_patch('epsilon1', 'epsilon1', limits, ax, b1=1)
    limit_patch('epsilon1', 'epsilon2', limits, ax, b1=3)
    limit_patch('epsilon2', 'epsilon1', limits, ax, b1=3)
    # limit_patch('epsilon2', 'epsilon2', limits, ax, b1=True)
    ax.plot(df[x][:-1], df[y][1:],  symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)
    ax.tick_params(labelleft='off')

    ax = ax_array[2, 2]
    x = 'zeta'
    y = 'delta'
    limit_patch('zeta1', 'delta1', limits, ax, b1=True)
    limit_patch('zeta2', 'delta2', limits, ax, b1=True)
    ax.plot(df[x], df[y],  symbol, c=gwp.qual_color(i_color))
    ax.tick_params(labelleft='off')  # , labelright='on')
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)

    return ax_array


if __name__ == '__main__':
    main()

    # Run pyclean
    try:
        subprocess.Popen(
            'pyclean .', shell=True, stderr=open('/dev/null', 'w'))
    except Exception:
        pass

    print 'completed get_dihedral script\n~~~~~ \m/ >.< \m/ ~~~~~~'




# dcd_files.append('raw_min.dcd')
# dcd_files.append('raw_min.dcd')
# dcd_files.append(run_dir + 'run0_1k_steps/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run1_1k_steps/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run2_100k_steps/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run2_100k_steps/monte_carlo/5p5_spb.dcd')
# dcd_files.append(run_dir + 'run3_100k_ngb/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run4_100k_ngb/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run5_100k_ngb/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run6_100k_ngb/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run7_100k_ngb/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run9_100k_ngb/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run10_100k_ngb/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run8_100k_ngb/monte_carlo/min_dsDNA60.dcd')
# dcd_files.append(run_dir + 'run3_100k_ngb/monte_carlo/min_dsDNA60_sparser.dcd')

# dcd_files.append(run_dir + 'run12_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run11_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run3_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run4_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run5_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run6_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run7_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run8_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run9_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
# dcd_files.append(run_dir + 'run10_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')

# dcd_files.append(run_dir + 'run12_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir + 'run11_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir + 'run3_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir + 'run4_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir + 'run5_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir + 'run6_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir + 'run7_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir + 'run8_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir +  [61, 120]]
# run_dir = '/home/schowell/data/myData/dihedrals/dsDNA_60bps/'
# pdb_file_name = run_dir + 'new_dsDNA60.pdb'
# flex_file = run_dir + 'new_dsDNA60.flex'
# drude = False
# dcd_files = []'run9_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')
# dcd_files.append(run_dir + 'run10_100k_ngb/energy_minimization/dsDNA_sparse_min.dcd')

# dcd_files.append(run_dir + 'run12_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run11_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run3_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run4_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run5_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run6_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run7_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run8_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run9_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run10_100k_ngb/energy_minimization_noMD/dsDNA_sparse_min_noMD.dcd')

# dcd_files.append(run_dir + 'run12_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run11_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run3_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run4_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run5_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run6_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run7_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run8_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run9_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run10_100k_ngb/c36_energy_minimization_1k_noMD/dsDNA_sparse_min_noMD.dcd')

# dcd_files.append(run_dir + 'run12_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run11_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run3_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run4_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run5_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run6_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run7_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run8_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run9_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
# dcd_files.append(run_dir + 'run10_100k_ngb/c36_energy_minimization_2k_noMD/dsDNA_sparse_min_noMD.dcd')
