#!/usr/bin/python
# Author:  Steven C. Howell
# Purpose: script for calculating DNA dihedral angles and generating plots
# Created: 22 September 2014

import sys
import os
import logging
import subprocess
import matplotlib
import multiprocessing as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import sassie.calculate.dna_dihedral as dd
import x_dna.util.gw_plot as gwp

LOGGER = logging.getLogger(__name__)  # add module name manually


class MainError(Exception):
    pass


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
    frequency = -1

    # ~~~~~~~~ FILE INPUT ~~~~~~~~~~ #
    # first and last DNA resdue of one of the DNA chains
    first_last_resids = [[1, 60], [61, 120]]
    run_dir = 'dsDNA_60bps/'
    pdb_file_name = run_dir + 'new_dsDNA60.pdb'
    flex_file = run_dir + 'new_dsDNA60.flex'

    processes = []
    scale_spb = 100.0/58.0
    max_frame = None
    used_goback = False
    show = False

    degree= u'\N{DEGREE SIGN}'
    # max_angles = np.array([1, 5, 10, 20, 30, 40, 50, 60])
    max_angles = np.array([60, 50, 40, 30, 20, 10, 5, 1])
    key = [r'$\delta\vartheta_{max}=%d$' % n + degree for n in max_angles]


    prefix = 'raw_MD'
    raw_dcd_files = []
    raw_dcd_files.append(run_dir + 'ma60_spb1000/mc1_spb1000_sparse.dcd')
    raw_dcd_files.append(run_dir + 'ma50_spb1000/mc1_spb1000_sparse.dcd')
    raw_dcd_files.append(run_dir + 'ma40_spb1000/mc1_spb1000_sparse.dcd')
    raw_dcd_files.append(run_dir + 'ma30_spb1000/mc1_spb1000_sparse.dcd')
    # raw_dcd_files.append(run_dir + 'ma20_spb1000/mc1_spb1000_sparse.dcd')
    raw_dcd_files.append(run_dir + 'ma20_spb1000/mc2_spb1000_sparse.dcd')
    raw_dcd_files.append(run_dir + 'ma10_spb1000/mc1_spb1000_sparse.dcd')
    raw_dcd_files.append(run_dir + 'ma5_spb1000/mc1_spb1000_sparse.dcd')
    raw_dcd_files.append(run_dir + 'ma1_spb1000/mc1_spb1000_sparse.dcd')
    inputs = (raw_dcd_files, key, max_frame, scale_spb, used_goback, show,
              prefix)
    processes.append(mp.Process(target=plot_all_dihedrals, args=inputs))
    dcd_files = raw_dcd_files

    prefix = 'mmm_MD'
    mmm_dcd_files = []
    mmm_dcd_files.append(run_dir + 'ma60_spb1000/mc1_spb1000_sparse_mmm.dcd')
    mmm_dcd_files.append(run_dir + 'ma50_spb1000/mc1_spb1000_sparse_mmm.dcd')
    mmm_dcd_files.append(run_dir + 'ma40_spb1000/mc1_spb1000_sparse_mmm.dcd')
    mmm_dcd_files.append(run_dir + 'ma30_spb1000/mc1_spb1000_sparse_mmm.dcd')
    # mmm_dcd_files.append(run_dir + 'ma20_spb1000/mc1_spb1000_sparse_mmm.dcd')
    mmm_dcd_files.append(run_dir + 'ma20_spb1000/mc2_spb1000_sparse_mmm.dcd')
    mmm_dcd_files.append(run_dir + 'ma10_spb1000/mc1_spb1000_sparse_mmm.dcd')
    mmm_dcd_files.append(run_dir + 'ma5_spb1000/mc1_spb1000_sparse_mmm.dcd')
    mmm_dcd_files.append(run_dir + 'ma1_spb1000/mc1_spb1000_sparse_mmm.dcd')
    inputs = (mmm_dcd_files, key, max_frame, scale_spb, used_goback, show,
              prefix)
    processes.append(mp.Process(target=plot_all_dihedrals, args=inputs))

    prefix = 'm_MD'
    m_dcd_files = []
    m_dcd_files.append(run_dir + 'ma60_spb1000/mc1_spb1000_sparse_m.dcd')
    m_dcd_files.append(run_dir + 'ma50_spb1000/mc1_spb1000_sparse_m.dcd')
    m_dcd_files.append(run_dir + 'ma40_spb1000/mc1_spb1000_sparse_m.dcd')
    m_dcd_files.append(run_dir + 'ma30_spb1000/mc1_spb1000_sparse_m.dcd')
    # m_dcd_files.append(run_dir + 'ma20_spb1000/mc1_spb1000_sparse_m.dcd')
    m_dcd_files.append(run_dir + 'ma20_spb1000/mc2_spb1000_sparse_m.dcd')
    m_dcd_files.append(run_dir + 'ma10_spb1000/mc1_spb1000_sparse_m.dcd')
    m_dcd_files.append(run_dir + 'ma5_spb1000/mc1_spb1000_sparse_m.dcd')
    m_dcd_files.append(run_dir + 'ma1_spb1000/mc1_spb1000_sparse_m.dcd')
    inputs = (m_dcd_files, key, max_frame, scale_spb, used_goback, show,
              prefix)
    processes.append(mp.Process(target=plot_all_dihedrals, args=inputs))

    # degree= u'\N{DEGREE SIGN}'
    # key = ['max: ' + str(n) + degree for n in np.arange(10, 50, 10)]
    # dcd_files = []
    # dcd_files.append(run_dir + 'run3_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
    # dcd_files.append(run_dir + 'run5_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
    # dcd_files.append(run_dir + 'run7_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
    # dcd_files.append(run_dir + 'run9_100k_ngb/monte_carlo/min_dsDNA60_sparse.dcd')
    # inputs = (dcd_files, key, max_frame, scale_spb, used_goback, show, prefix)
    # processes.append(mp.Process(target=plot_all_dihedrals, args=inputs))
    # ~~~~~~~~ END FILE INPUT ~~~~~~~~~~ #

    # ~~~~~~~~ RUN INPUT ~~~~~~~~~~ #

    # ~~~~~~~~ RUN INPUT ~~~~~~~~~~ #

    # plot_aez_dihedrals(dcd_files, key, max_frame=None, scale_spb=(127.0/58.0),
                        # show=False, prefix=prefix)
    #
    # processes = []
    # if processes:
        # for p in processes:
            # p.start()

        # for p in processes:
            # p.join()
    # else:
        # plot_all_dihedrals(dcd_files, key, max_frame=None,
                           # scale_spb=(127.0/58.0), show=False, prefix=prefix)

    compare_min_dihedrals(raw_dcd_files, m_dcd_files, key, show=False,
                          scale_spb=scale_spb, prefix='c36_MD_m')
    compare_min_dihedrals(raw_dcd_files, mmm_dcd_files, key, show=False,
                          scale_spb=scale_spb, prefix='c36_MD_mmm')


def compare_min_dihedrals(raw_dcd_files, min_dcd_files, key, max_frame=None,
                       scale_spb=1, used_goback=False, show=False, prefix=None):
    raw_angle_vs_steps = []
    min_angle_vs_steps = []
    raw_all_in_range = []
    min_all_in_range = []

    for dcd_file_name in min_dcd_files:
        if max_frame:
            npy_file = dcd_file_name[:-4] + '_%d.npy' % max_frame
        else:
            npy_file = dcd_file_name[:-3] + 'npy'

        all_in_range_file = dcd_file_name[:-4] + '_ave_in_range'
        try:
            if 'MD' in prefix:
                npy_file = npy_file.replace('.npy', '_MD.npy')
                all_in_range_file += '_MD'
        except:
            pass

        min_angle_vs_steps.append(np.load(npy_file))
        min_all_in_range.append(np.loadtxt(all_in_range_file + '.out'))
        print 'loaded percent in range files: \n%s \n%s' % (npy_file,
                                                            all_in_range_file)

    for dcd_file_name in raw_dcd_files:
        if max_frame:
            npy_file = dcd_file_name[:-4] + '_%d.npy' % max_frame
        else:
            npy_file = dcd_file_name[:-3] + 'npy'

        all_in_range_file = dcd_file_name[:-4] + '_ave_in_range'
        try:
            if 'MD' in prefix:
                npy_file = npy_file.replace('.npy', '_MD.npy')
                all_in_range_file += '_MD'
        except:
            pass

        raw_angle_vs_steps.append(np.load(npy_file))
        raw_all_in_range.append(np.loadtxt(all_in_range_file + '.out'))
        print 'loaded percent in range files: \n%s \n%s' % (npy_file,
                                                            all_in_range_file)

    fig = plt.figure(figsize=(5, 12))
    gs = gridspec.GridSpec(8, 2, left=0.1, right=0.75, wspace=0.07, hspace=0)
    angle_labels = [r'$\alpha$', r'$\beta$', r'$\gamma$', r'$\delta$',
                    r'$\epsilon$', r'$\zeta$', r'$\chi$', r'All']
    plot_order = [0, 2, 3, 4, 7, 1, 5, 6] # change to show select angles
    default_fontsize = 12

    def plot_settings(ax, x_min, x_max):
        ax.set_ylim([71, 104])
        ax.set_xscale('log')
        ax.set_xlim([min(x_min), min(x_max)])

    for (i, i_angle) in enumerate(plot_order):
        ax = plt.subplot(gs[i, 0])
        plt.axhline(95, ls='--', c='DimGray', linewidth=2)
        x_max = []
        x_min = []
        for j in xrange(len(raw_dcd_files)):
            if i_angle == 0:
                ax.plot(raw_all_in_range[j][:, 0],
                         raw_all_in_range[j][:, 1]*100, '-',
                         c=gwp.qual_color(j), label=key[j])
                x_max.append(raw_all_in_range[j][-1, 0])
                x_min.append(raw_all_in_range[j][0, 0])
                ax.text(18, 107, 'Raw')
            else:
                ax.plot(raw_angle_vs_steps[j][:, 0],
                         raw_angle_vs_steps[j][:, i_angle]*100, '-',
                         c=gwp.qual_color(j), label=key[j])
                x_max.append(raw_angle_vs_steps[j][-1, 0])
                x_min.append(raw_angle_vs_steps[j][0, 0])
        plot_settings(ax, x_min, x_max)
        ax.set_ylabel(r'% in range')
        if i < len(plot_order) - 1:
            ax.xaxis.set_ticklabels([])
            # ax1.xaxis.set_ticks([])
        else:
            ax.set_xlabel(r'steps per bp')
        ax.set_title(angle_labels[i_angle-1], y=0.03, x=0.12)

    for (i, i_angle) in enumerate(plot_order):
        ax = plt.subplot(gs[i, 1])
        plt.axhline(95, ls='--', c='DimGray', linewidth=2)
        for j in xrange(len(min_dcd_files)):
            if i_angle == 0:
                ax.plot(min_all_in_range[j][:, 0],
                         min_all_in_range[j][:, 1]*100, '-',
                         c=gwp.qual_color(j), label=key[j])
                ax.text(9, 107, 'Minimized')
            else:
                ax.plot(min_angle_vs_steps[j][:, 0],
                         min_angle_vs_steps[j][:, i_angle]*100, '-',
                         c=gwp.qual_color(j), label=key[j])
        plot_settings(ax, x_min, x_max)
        ax.yaxis.set_ticklabels([])
        if i < len(plot_order) - 1:
            ax.xaxis.set_ticklabels([])
            # ax1.xaxis.set_ticks([])
        else:
            ax.set_xlabel(r'steps per bp')
        ax.set_title(angle_labels[i_angle-1], y=0.03, x=0.12)

        if i == 0:
            lg = plt.legend(loc='upper left', bbox_to_anchor=(1, 1),
                            scatterpoints=1, numpoints=1,
                            prop={'size': default_fontsize})

    path = 'dsDNA_60bps/'
    save_name = 'raw_min_percent_in_range'
    if prefix:
        save_name = prefix + '_' + save_name
    if max_frame:
        save_name += '_' + str(max_frame)
    plt.savefig(path + save_name + '.eps', dpi=400, bbox_inches='tight')
    plt.savefig(path + save_name + '.png', dpi=400, bbox_inches='tight')
    print 'saved figures to %s.eps and %s.png' % (save_name, save_name)
    if show:
        plt.show()

    print 'done plotting'

def plot_all_dihedrals(dcd_files, key, max_frame=None, scale_spb=1,
                       used_goback=False, show=False, prefix=None):
    angle_vs_steps = []
    for dcd_file_name in dcd_files:
        if max_frame:
            npy_file = dcd_file_name[:-4] + '_%d.npy' % max_frame
        else:
            npy_file = dcd_file_name[:-3] + 'npy'

        try:
            if 'MD' in prefix:
                npy_file = npy_file.replace('.npy', '_MD.npy')
        except:
            pass

        angle_vs_steps.append(np.load(npy_file))
        print 'loaded percent of good angles file: %s' % npy_file


    fig = plt.figure(figsize=(4, 20))
    gs1 = gridspec.GridSpec(7, 1, left=0.1, right=0.75, wspace=0, hspace=0)
    angle_labels = [r'$\alpha$', r'$\beta$', r'$\gamma$', r'$\delta$',
                    r'$\epsilon$', r'$\zeta$', r'$\chi$']
    plot_order = [2, 3, 4, 7, 1, 5, 6]  # only plot alpha, beta, and gamma
    default_fontsize = 12
    for (i, i_angle) in enumerate(plot_order):
        ax1 = plt.subplot(gs1[i])
        x_max = []
        for j in xrange(len(dcd_files)):
            ax1.plot(angle_vs_steps[j][:, 0],
                     angle_vs_steps[j][:, i_angle]*100, '-',
                     c=gwp.qual_color(j), label=key[j])
            x_max.append(angle_vs_steps[j][-1, 0])

        ax1.set_ylim([73, 102])
        ax1.set_ylabel(r'% in range')
        ax1.set_xscale('log')
        ax1.set_xlim([0, min(x_max)])
        if i < len(plot_order) - 1:
            ax1.xaxis.set_ticklabels([])
            # ax1.xaxis.set_ticks([])
        else:
            ax1.set_xlabel(r'steps per bp')
        if i == 0:
            lg = plt.legend(loc='upper left', bbox_to_anchor=(1, 1),
                            scatterpoints=1, numpoints=1,
                            prop={'size': default_fontsize})
        ax1.set_title(angle_labels[i_angle-1], y=0.05)

    if show:
        plt.show()
    else:
        save_name = 'dsDNA_60bps/all_in_range'
        if prefix:
            save_name = save_name + '_' + prefix
        if max_frame:
            save_name += '_' + str(max_frame)
        plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
        plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
        print 'saved figures to %s.eps and %s.png' % (save_name, save_name)
    print 'done plotting'


def get_step_i_dihedrals(i, indices, all_df):
    df_list = []
    frames = np.where(indices == i)[0]
    for i_angle in frames:
        df_list.append(all_df[all_df['frame'] == i_angle + 1])
    i_df = pd.concat(df_list)
    return i_df, len(frames)



def plot_aez_dihedrals(dcd_files, key, max_frame=None, scale_spb=1,
                       used_goback=False, show=False, prefix=None):
    angle_vs_steps = []
    for dcd_file_name in dcd_files:
        if max_frame:
            npy_file = dcd_file_name[:-4] + '_%d.npy' % max_frame
        else:
            npy_file = dcd_file_name[:-3] + 'npy'

        angle_vs_steps.append(np.load(npy_file))
        print 'loaded percent of good angles file: %s' % npy_file


    fig = plt.figure(figsize=(6.5, 10))
    gs1 = gridspec.GridSpec(3, 1, left=0.1, right=0.75, wspace=0, hspace=0)
    angle_labels = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi']
    plot_order = [1, 5, 6]  # only plot alpha, beta, and gamma
    default_fontsize = 12
    for (i, i_angle) in enumerate(plot_order):
        ax1 = plt.subplot(gs1[i])
        x_max = []
        for j in xrange(len(dcd_files)):
            ax1.plot(angle_vs_steps[j][:, 0], angle_vs_steps[j][:, i_angle], '-',
                     c=gwp.qual_color(j), label=key[j])
            x_max.append(angle_vs_steps[j][-1, 0])

        ax1.set_ylim([0, 1.1])
        ax1.set_ylabel(r'% in range')
        ax1.set_xlim([0, min(x_max)])
        if i < len(plot_order) - 1:
            ax1.xaxis.set_ticklabels([])
            # ax1.xaxis.set_ticks([])
        else:
            ax1.set_xlabel(r'steps per bp')
        if i == 0:
            lg = plt.legend(loc='upper left', bbox_to_anchor=(1, 1),
                            scatterpoints=1, numpoints=1,
                            prop={'size': default_fontsize})
        ax1.set_title(angle_labels[i_angle-1], y=0.05)

    if show:
        plt.show()
    else:
        save_name = 'dsDNA_60bps/moved_in_range'
        if prefix:
            save_name = save_name + '_' + prefix
        if max_frame:
            save_name += '_' + str(max_frame)
        plt.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
        plt.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
        print 'saved figures to %s.eps and %s.png' % (save_name, save_name)
    print 'done plotting'


def limit_patch(x_key, y_key, limits, ax):
    red_alpha_p5 = np.array([235, 157, 158], dtype=float) / 255
    limit_patch = patches.Rectangle(
        (limits.loc[x_key]['low'], limits.loc[y_key]['low']),
        4 * limits.loc[x_key]['sd'], 4 * limits.loc[y_key]['sd'],
        edgecolor='none', alpha=1, facecolor=red_alpha_p5)
    ax.add_patch(limit_patch)


def scatter_plot_dihedrals(dcd_file_name, frequency=-1, scale_spb=1):
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

    limits = pd.DataFrame.from_csv('~/data/myData/dihedrals/'
                                   'dsDNA_60bps/dna_angles.txt')
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
        df_0, _ = get_step_i_dihedrals(1, indices, angles)
        ax_array = make_selected_scatter_plots(limits, df_0, angle_label,
                                               symbol='+')

        df_1, _ = get_step_i_dihedrals(0, indices, angles)
        make_selected_scatter_plots(limits, df_1, angle_label, i_color=2,
                                    ax_array=ax_array)

        plt.legend(['raw', 'min'], loc='upper left',
                   bbox_to_anchor=[1.1475, 1.07], numpoints=1)

        plt.suptitle('Comparative scatter plot of selected torsional angles')
        plt.savefig(prefix + 'dihedrals_comparison.eps', dpi=400,
                    bbox_inches='tight')
        plt.savefig(prefix + 'dihedrals_comparison.png', dpi=400,
                    bbox_inches='tight')
        # plt.show()
    else:
        if frequency == -1:
            frequency = n_steps / 10
        for step in xrange(frequency-1, n_steps, frequency):
            df, nf = get_step_i_dihedrals(step, indices, angles)
            make_selected_scatter_plots(limits, df, angle_label)

            plt.suptitle('%d Steps (%d frame/s): scatter plots of selected'
                         'torsional angles' % (step + 1, nf))
            plt.savefig(prefix + 'dihedrals_%dsteps.eps' %
                        (step + 1), dpi=400, bbox_inches='tight')
            plt.savefig(prefix + 'dihedrals_%dsteps.png' %
                        (step + 1), dpi=400, bbox_inches='tight')

    # plt.show()
    print 'pause'

def plt_format(ax):

    ax.set_xlim([0, 360])
    ax.set_ylim([0, 360])
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_yticks([0, 60, 120, 180, 240, 300, 360])

def make_selected_scatter_plots(limits, df, angle_label, ax_array=[],
                                symbol='x', i_color=0):

    if len(ax_array) == 0:
        fig, ax_array = plt.subplots(3, 3)

    ax = ax_array[0, 0]
    x = 'zeta'
    y = 'alpha'
    limit_patch('zeta1', 'alpha', limits, ax)
    limit_patch('zeta2', 'alpha', limits, ax)
    ax.plot(df[x][0:-1], df[y][1:], symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)
    ax.tick_params(labelbottom='off') # , labeltop='on')

    ax = ax_array[0, 1]
    x = 'zeta'
    y = 'beta'
    limit_patch('zeta1', 'beta1', limits, ax)
    limit_patch('zeta1', 'beta2', limits, ax)
    limit_patch('zeta2', 'beta1', limits, ax)
    limit_patch('zeta2', 'beta2', limits, ax)
    ax.plot(df[x][0:-1], df[y][1:], symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)
    ax.tick_params(labelleft='off', labelbottom='off') # , labeltop='on')

    ax = ax_array[0, 2]
    x = 'zeta'
    y = 'epsilon'
    limit_patch('zeta1', 'epsilon1', limits, ax)
    limit_patch('zeta2', 'epsilon2', limits, ax)
    ax.plot(df[x], df[y], symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    ax.tick_params(
        labelleft='off', labelbottom='off')# ,  labeltop='on', labelright='on')

    ax = ax_array[1, 0]
    x = 'alpha'
    y = 'gamma'
    limit_patch('alpha', 'gamma', limits, ax)
    ax.plot(df[x], df[y], symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    ax.tick_params(labelbottom='off')

    ax = ax_array[1, 1]
    x = 'zeta'
    y = 'chi'
    limit_patch('zeta1', 'chi1pu', limits, ax)
    limit_patch('zeta1', 'chi1py', limits, ax)
    limit_patch('zeta2', 'chi2', limits, ax)
    ax.plot(df[x], df[y],  symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    ax.tick_params(labelleft='off', labelbottom='off')

    ax = ax_array[1, 2]
    x = 'delta'
    y = 'chi'
    limit_patch('delta1', 'chi1pu', limits, ax)
    limit_patch('delta1', 'chi1py', limits, ax)
    limit_patch('delta2', 'chi2', limits, ax)
    ax.plot(df[x], df[y],  symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y])
    plt_format(ax)
    # , labelright='on')
    ax.tick_params(labelleft='off', labelbottom='off')

    ax = ax_array[2, 0]
    x = 'zeta'
    y = 'zeta'
    limit_patch('zeta1', 'zeta1', limits, ax)
    limit_patch('zeta1', 'zeta2', limits, ax)
    limit_patch('zeta2', 'zeta1', limits, ax)
    # limit_patch('zeta2', 'zeta2', limits, ax)
    ax.plot(df[x][:-1], df[y][1:],  symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)

    ax = ax_array[2, 1]
    x = 'epsilon'
    y = 'epsilon'
    limit_patch('epsilon1', 'epsilon1', limits, ax)
    limit_patch('epsilon1', 'epsilon2', limits, ax)
    limit_patch('epsilon2', 'epsilon1', limits, ax)
    # limit_patch('epsilon2', 'epsilon2', limits, ax)
    ax.plot(df[x][:-1], df[y][1:],  symbol, c=gwp.qual_color(i_color))
    ax.set_xlabel(angle_label[x])
    ax.set_ylabel(angle_label[y] + ' + 1')
    plt_format(ax)
    ax.tick_params(labelleft='off')

    ax = ax_array[2, 2]
    x = 'zeta'
    y = 'delta'
    limit_patch('zeta1', 'delta1', limits, ax)
    limit_patch('zeta2', 'delta2', limits, ax)
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
