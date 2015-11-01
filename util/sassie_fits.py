#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Compare experimental data to the sassie structures
# Created: 20 March 2015
#
# $Id$
#
# 000000001111111111222222222233333333334444444444555555555566666666667777777777
# 234567890123456789012345678901234567890123456789012345678901234567890123456789


import logging
LOGGER = logging.getLogger(__name__)  # add module name manually

import os
import glob
import locale
import errno
import shutil
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy import optimize
try:
    dummy = os.environ["DISPLAY"]
except:
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import sassie.sasmol.sasmol as sasmol
import x_dna.util.gw_plot as gp
# import x_dna.drivers.myAlign as align
import x_dna.util.plt_inset as plt_inset
debug = True
op = os.path


class MainError(Exception):
    pass


def pr():
    NotImplemented


def load_crysol(saspath, i0):
    '''
    load in the crysol output (Rg and I(Q))
    taken from "sassie/analyze/best.py"
    '''

    sasfilestring = saspath.split('/')[-2] + '_'

    nfilestring = 'find ' + saspath + ' -name "*.log" | grep -c log'
    nfout = os.popen(nfilestring, 'r').readlines()
    nf = locale.atoi(nfout[0])

    print "\n# READING Rg VALUES FROM LOG FILES"
    log = []
    assert nf > 0, 'no I(Q) calculations in %s' % saspath
    for i in xrange(nf):

        mst = str(i + 1).zfill(5)  # 99999 maximum number of frames
        log.append(saspath + '/' + sasfilestring + mst)

    Rg = crysol_rg(nf, log)
    iq_data = load_crysol_iq(nf, log, i0)

    return Rg, iq_data, log


def load_foxs(saspath):
    '''
    load in the FoXS output (Rg and I(Q))
    '''
    result_file = op.join(saspath, 'rg.csv')

    syn_files = glob.glob(op.join(saspath, '*.dat'))
    assert len(syn_files) > 0, 'no I(Q) calculations in %s' % saspath
    syn_files.sort()

    rg = [None] * len(syn_files)
    if op.isfile(result_file):
        rg_df = pd.read_csv(result_file, sep='\t')
        rg_df.index = rg_df['labels']
        save_rg = False
        for (i, syn_file) in enumerate(syn_files):
            rg[i] = rg_df['rg'].loc[op.split(syn_file)[1].split('.')[0]]
    else:
        try:
            dcd_file = glob.glob(op.join(op.join(op.split(saspath)[0],
                                                 'monte_carlo'), '*.dcd'))
            assert len(dcd_file) == 1, 'ERROR: not clear which dcd file to use'
            dcd_file = dcd_file[0]
            pdb_search = op.join(op.split(op.split(saspath)[0])[0], '*.pdb')
            pdb_file = glob.glob(pdb_search)[0]
            print '\ncalculating the Rg from this dcd: %s' % dcd_file
            print 'together with this pdb file: %s' % pdb_file
            rg = dcd_rg(pdb_file, dcd_file)
        except:
            print 'did not calculate Rg using dcd and pdb'
        save_rg = True

    # all_data = []
    iq = []
    q = []
    label = []
    for (i, syn_file) in enumerate(syn_files):
        data = np.loadtxt(syn_file)
        # all_data.append(data)
        iq.append(data[:, 1])
        q.append(data[:, 0])
        label.append(op.basename(syn_file).split('.')[0])
        if not rg[i]:
            pdb_file = syn_file.replace('foxs', 'pdb').replace('.dat', '')
            rg[i] = pdb_rg(pdb_file)

    if save_rg:
        rg_dict = {'rg': rg}
        rg_df = pd.DataFrame(rg_dict, index=label)
        rg_df.index.name = 'labels'
        rg_df.to_csv(result_file, sep='\t')
    rg = np.array(rg)

    iq_data = np.concatenate(
        (q[0][..., None], np.array(iq).transpose()), axis=1)

    return rg, iq_data, label


def dcd_rg(pdb_file, dcd_file):
    '''
    given a pdb filename, this will return the radius of gyration
    '''
    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_file)
    mol.read_dcd(dcd_file)
    rg = []
    for i in xrange(mol.number_of_frames()):
        rg.append(mol.calcrg(i))
    return rg


def pdb_rg(pdb_file):
    '''
    given a pdb filename, this will return the radius of gyration
    '''
    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_file)
    return mol.calcrg(0)


def load_crysol_iq(ns, log, i0=None):
    '''
    read in the I(q) from the crysol output
    taken from "sassie/analyze/best.py"
    '''
    norm = 1.0
    for i in xrange(ns):
        sinf = open(log[i] + '.int', 'r').readlines()
        nh = 1
        sln = len(sinf)
        jj = 0
        if 0 == i:
            # on the first iteration, setup the output array
            allspec = np.zeros((ns, (sln - nh)), np.float)
            qvar = np.zeros((sln - nh), np.float)
            print '# READING ' + str(ns) + ' SAS INT FILES'
        for j in xrange(sln):
            if i0:
                slin = sinf[j].split()
                if(j >= nh):
                    qval = locale.atof(slin[0])
                    qvar[jj] = qval
                    val1 = locale.atof(slin[1])
                    if(jj == 0):
                        norm = val1
                    nval1 = i0 * val1 / norm
                    allspec[i][jj] = nval1
                    jj = jj + 1
            else:
                slin = sinf[j].split()
                if(j >= nh):
                    qvar[jj] = locale.atof(slin[0])
                    allspec[i][jj] = locale.atof(slin[1])
                    jj += 1
    # not sure pandas is the way to go
    # iq_data = DataFrame(allspec, columns=qvar, index=range(1,ns+1)).T
    # iq_data.columns.name = 'id' ; iq_data.index.name = 'Q'

    # numpy seem to be the best option the for storing the I(Q) data
    iq_data = np.zeros(((sln - nh), ns + 1))
    iq_data[:, 0] = qvar
    iq_data[:, 1:] = allspec.T

    return iq_data


def gnom_rg(gnom_files):
    '''
    read the Rg in from the GNOM output
    modified from "sassie/analyze/best.py"
    '''
    prefix = []
    rg_reci = []
    rg_real = []
    rger_real = []
    i0_reci = []
    i0_real = []
    i0er_real = []
    for gnom_out in gnom_files:
        inf = open(gnom_out, 'r').readlines()
        ln = len(inf)
        for k in xrange(ln):
            lin = inf[k].split()
            # print lin
            if(len(lin) > 0):
                if(lin[0] == 'Reciprocal'):
                    rg_reci.append(locale.atof(lin[4], func=float))
                    i0_reci.append(locale.atof(lin[8], func=float))
                elif(lin[0] == 'Real' and lin[1] == 'space:'):
                    rg_real.append(locale.atof(lin[4], func=float))
                    rger_real.append(locale.atof(lin[6], func=float))
                    i0_real.append(locale.atof(lin[9], func=float))
                    i0er_real.append(locale.atof(lin[11], func=float))
                    index = op.split(op.splitext(gnom_out)[0])[-1]
                    prefix.append(index.replace('_zeroCon', ''))

    gnom_dict = {'Rg gq': rg_reci, 'Rg gr': rg_real, 'RgEr gr': rger_real,
                 'I0 gq': i0_reci, 'I0 gr': i0_real, 'I0Er gr': i0er_real}
    return pd.DataFrame(gnom_dict, index=prefix)


def crysol_rg(ns, log):
    '''
    read the Rg in from the crysol output
    taken from "sassie/analyze/best.py"
    '''
    rgarray = []
    keep = []
    rgcl = []
    rgch = []
    for i in xrange(ns):
        inf = open(log[i] + '.log', 'r').readlines()
        ln = len(inf)
        for k in xrange(ln):
            lin = inf[k].split()
            if(len(lin) > 0):
                if(lin[0] == 'Electron'):
                    trg = lin[3]
                    ftrg = locale.atof(trg)
                    keep.append(i)
                    rgarray.append(ftrg)

    # rg_series = Series(data=rgarray, index=range(1,ns+1))
    # return rg_series
    return np.array(rgarray)

# def load_iq(saspath, ns):
    # '''
    # read in the crysol I(Q) output
    # taken from "sassie/analyze/best.py"
    # '''


def mkdir_p(path):
    '''
    make directory recursively
    adapted from http://stackoverflow.com/questions/600268/
    '''
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and op.isdir(path):
            pass
        else:
            raise


def compare_run_to_iq_all_metrics(run_dir, goal_iq, filter_dir, out_file=None,
                           run_label='', match_range=None):
    '''
    compare the data and models using all discrepancy metrics
    '''
    assert op.exists(run_dir), 'No such run directory: %s' % run_dir
    if op.exists(op.join(run_dir, 'crysol')) and False:
        syn_data_dir = op.join(run_dir, 'crysol')
        ext = '/*.int'
    elif op.exists(op.join(run_dir, 'foxs')):
        syn_data_dir = op.join(run_dir, 'foxs')
        ext = '/*.dat'
    else:
        assert False, 'failed to find the calculated scattering data'

    # get the Rg and I(Q) for each structure from the calculation rusults
    if ext[-3:] == 'int':
        Rg, data_iq, labels = load_crysol(syn_data_dir, goal_iq[0, 1])
    elif ext[-3:] == 'dat':
        Rg, data_iq, labels = load_foxs(syn_data_dir) #, goal_iq[0, 1])

    data_iq = interp_iq(data_iq, goal_iq)

    # compare I(Q) of the experimental and synthetic structures
    nf = data_iq.shape[1] - 1
    x2_scaled_iq = np.zeros(data_iq.shape)
    x2_scaled_iq[:, 0] = data_iq[:, 0]
    wr_scaled_iq = np.copy(x2_scaled_iq)
    i1_scaled_iq = np.copy(x2_scaled_iq)
    mx2_scale = np.zeros(nf)
    mwr_scale = np.zeros(nf)
    i1_scale = np.zeros(nf)
    iq_offset = np.zeros(nf)

    # X^2
    mx2_x2 = np.zeros(nf) # minimizing the x2
    mwr_x2 = np.zeros(nf) # minimizing the weighted R-factor
    i1_x2 = np.zeros(nf)  # scaling to the first data point

    # weighted R-factor
    mx2_wr = np.zeros(nf) # minimizing the x2
    mwr_wr = np.zeros(nf) # minimizing the weighted R-factor
    i1_wr = np.zeros(nf)  # scaling to the first data point

    # R-factor
    mx2_r = np.zeros(nf) # minimizing the x2
    mwr_r = np.zeros(nf) # minimizing the weighted R-factor
    i1_r = np.zeros(nf)  # scaling to the first data point

    # F-factor
    mx2_f = np.zeros(nf) # minimizing the x2
    mwr_f = np.zeros(nf) # minimizing the weighted R-factor
    i1_f = np.zeros(nf)  # scaling to the first data point

    for i in xrange(1, nf+1):
        j = i - 1
        # if s and o:
            # match_iq[:, [0, i]], iq_scale[j], iq_offset[j], mx2_x2[j] = scale_offset(
                # data_iq[:, [0, i]], goal_iq)

        iq_offset[j] = 0
        if match_range:
            # only scale a specific Q-range (e.g. high-Q range)
            i_min = np.argmin(np.abs(goal_iq[:, 0] - match_range[0]))
            i_max = np.argmin(np.abs(goal_iq[:, 0] - match_range[1]))

            # minimize chi^2
            _, mx2_scale[j], _ = scale(data_iq[i_min:i_max+1, [0, i]],
                                       goal_iq[i_min:i_max+1, :])
            x2_scaled_iq[:, i] = mx2_scale[j] * data_iq[:, i]
            mx2_x2[j] = get_x2(goal_iq, x2_scaled_iq[:, [0, i]])

            # minimize weighted R-factor
            _, mwr_scale[j], _ = scale_wr(data_iq[i_min:i_max+1, [0, i]],
                                          goal_iq[i_min:i_max+1, :])
            wr_scaled_iq[:, i] = mwr_scale[j] * data_iq[:, i]
            mwr_wr[j] = get_wr(goal_iq, wr_scaled_iq[:, [0, i]])

        else:
            # match the entire Q-range
            x2_scaled_iq[:, [0, i]], mx2_scale[j], mx2_x2[j] = scale(
                data_iq[:, [0, i]], goal_iq)
            wr_scaled_iq[:, [0, i]], mwr_scale[j], mwr_wr[j] = scale_wr(
                data_iq[:, [0, i]], goal_iq)

        # scale using the first data point
        i1_scale[j] = goal_iq[0, 1]/data_iq[0, i]
        i1_scaled_iq[:, i] = i1_scale[j] * data_iq[:, i]

        # X^2
        mwr_x2[j] = get_x2(goal_iq, wr_scaled_iq[:, [0, i]])
        i1_x2[j] = get_x2(goal_iq, i1_scaled_iq[:, [0, i]])

        # weighted R-factor
        mx2_wr[j] = get_wr(goal_iq, x2_scaled_iq[:, [0, i]])
        i1_wr[j] = get_wr(goal_iq, i1_scaled_iq[:, [0, i]])

        # R-factor
        mx2_r[j] = get_r(goal_iq, x2_scaled_iq[:, [0, i]])
        mwr_r[j] = get_r(goal_iq, wr_scaled_iq[:, [0, i]])
        i1_r[j] = get_r(goal_iq, i1_scaled_iq[:, [0, i]])

        # F-factor
        mx2_f[j] = get_f(goal_iq, x2_scaled_iq[:, [0, i]])
        mwr_f[j] = get_f(goal_iq, wr_scaled_iq[:, [0, i]])
        i1_f[j] = get_f(goal_iq, i1_scaled_iq[:, [0, i]])

    res_dict = {'Rg': Rg, 'offset': iq_offset, 'labels': labels,
                'mx2_scale': mx2_scale, 'mwr_scale': mwr_scale, 'i1_scale': i1_scale,
                'mx2_x2': mx2_x2, 'mwr_x2': mwr_x2, 'i1_x2': i1_x2,
                'mx2_wr': mx2_wr, 'mwr_wr': mwr_wr, 'i1_wr': i1_wr,
                'mx2_r': mx2_r,   'mwr_r': mwr_r,   'i1_r': i1_r,
                'mx2_f': mx2_f,   'mwr_f': mwr_f,   'i1_f': i1_f}
    result_df = pd.DataFrame(res_dict, index=range(1, nf+1))
    result_df.index.name = 'id'

    # save output to filter directory
    mkdir_p(filter_dir)
    if not out_file:
        out_file = op.join(filter_dir, 'rg_x2.out')
    result_df.to_csv(out_file, float_format='%5.10f', sep='\t')

    # small file
    np.savetxt(op.join(filter_dir, 'goal%s.iq' % run_label), goal_iq)
    # too big to be useful as text file
    np.save(op.join(filter_dir, 'x2_data_iq%s.npy' % run_label), x2_scaled_iq)
    np.save(op.join(filter_dir, 'wr_data_iq%s.npy' % run_label), wr_scaled_iq)
    np.save(op.join(filter_dir, 'i1_data_iq%s.npy' % run_label), i1_scaled_iq)

    return result_df, x2_scaled_iq, wr_scaled_iq, i1_scaled_iq, goal_iq


def compare_run_to_iq(run_dir, goal_iq, filter_dir, metric='r', out_file=None,
                           run_label='', match_range=None):
    '''
    compare the data and models using selected discrepancy metrics
    '''
    assert op.exists(run_dir), 'No such run directory: %s' % run_dir
    if op.exists(op.join(run_dir, 'crysol')) and False:
        syn_data_dir = op.join(run_dir, 'crysol')
        ext = '/*.int'
    elif op.exists(op.join(run_dir, 'foxs')):
        syn_data_dir = op.join(run_dir, 'foxs')
        ext = '/*.dat'
    else:
        assert False, 'failed to find the calculated scattering data'

    # get the Rg and I(Q) for each structure from the calculation rusults
    if ext[-3:] == 'int':
        Rg, data_iq, labels = load_crysol(syn_data_dir, goal_iq[0, 1])
    elif ext[-3:] == 'dat':
        Rg, data_iq, labels = load_foxs(syn_data_dir) #, goal_iq[0, 1])

    data_iq = interp_iq(data_iq, goal_iq)

    # setup the output storage
    nf = data_iq.shape[1] - 1

    i1_scaled_iq = np.zeros(data_iq.shape)
    i1_scaled_iq[:, 0] = data_iq[:, 0]
    min_scaled_iq = np.copy(i1_scaled_iq)

    i1_scale = np.zeros(nf)
    min_scale = np.zeros(nf)

    iq_offset = np.zeros(nf)
    min_offset = np.zeros(nf)

    i1_residual = np.zeros(nf) # scaling to the first data point
    min_residual = np.zeros(nf) # minimizing the discrepancy

    for i in xrange(1, nf+1):
        j = i - 1
        # if s and o:
            # match_iq[:, [0, i]], iq_scale[j], iq_offset[j], mx2_x2[j] = scale_offset(
                # data_iq[:, [0, i]], goal_iq)

        iq_offset[j] = 0

        if metric == 'r':
                # do not minimize the discrepancy
                min_scale = i1_scale
                min_scaled_iq = i1_scaled_iq
                min_residual = i1_residual

        elif match_range:
            # only scale a specific Q-range (e.g. high-Q range)
            i_min = np.argmin(np.abs(goal_iq[:, 0] - match_range[0]))
            i_max = np.argmin(np.abs(goal_iq[:, 0] - match_range[1]))

            if metric == 'x2':
                # minimize chi^2
                _, min_scale[j], _ = scale(data_iq[i_min:i_max+1, [0, i]],
                                           goal_iq[i_min:i_max+1, :])
                min_scaled_iq[:, i] = min_scale[j] * data_iq[:, i]
                min_residual[j] = get_x2(goal_iq, min_scaled_iq[:, [0, i]])
            elif metric == 'wR':
                # minimize weighted R-factor
                _, min_scale[j], _ = scale_wr(data_iq[i_min:i_max+1, [0, i]],
                                              goal_iq[i_min:i_max+1, :])
                min_scaled_iq[:, i] = min_scale[j] * data_iq[:, i]
                min_residual[j] = get_wr(goal_iq, min_scaled_iq[:, [0, i]])
            else:
                assert False, ('ERROR: use predefined metric (x2, wr, or r), '
                               'or define a new one')
        else:
            # match the entire Q-range
            if metric == 'x2':
                min_scaled_iq[:, [0, i]], min_scale[j], min_residual[j] = scale(
                    data_iq[:, [0, i]], goal_iq)
            elif metric == 'wR':
                min_scaled_iq[:, [0, i]], min_scale[j], min_residual[j] = scale_wr(
                    data_iq[:, [0, i]], goal_iq)
            else:
                assert False, ('ERROR: use predefined metric (x2, wr, or r), '
                               'or define a new one')

        # scale using the first data point
        i1_scale[j] = goal_iq[0, 1]/data_iq[0, i]
        i1_scaled_iq[:, i] = i1_scale[j] * data_iq[:, i]

        if metric == 'x2':
            # X^2
            i1_residual[j] = get_x2(goal_iq, i1_scaled_iq[:, [0, i]])
            min_key = 'min_%s' % metric
        elif metric == 'wR':
            # weighted R-factor
            i1_residual[j] = get_wr(goal_iq, i1_scaled_iq[:, [0, i]])
            min_key = 'min_%s' % metric
        elif metric == 'r':
            # R-factor
            i1_residual[j] = get_r(goal_iq, i1_scaled_iq[:, [0, i]])
            min_key = 'min_i0'

    i1_key = 'i1_%s' % metric
    res_dict = {'Rg': Rg, 'offset': iq_offset, 'labels': labels, 'min_scale':
                min_scale, min_key: min_residual, i1_key: i1_residual}
    result_df = pd.DataFrame(res_dict, index=range(1, nf+1))
    result_df.index.name = 'id'

    # save output to filter directory
    mkdir_p(filter_dir)
    if not out_file:
        out_file = op.join(filter_dir, 'rg_residual.out')
    result_df.to_csv(out_file, float_format='%5.10f', sep='\t')

    # small file
    np.savetxt(op.join(filter_dir, 'goal_%s.iq' % run_label), goal_iq)

    # too big to be useful as text file
    np.save(op.join(filter_dir, 'min_data_%s.npy' % run_label), min_scaled_iq)
    np.save(op.join(filter_dir, 'i1_data_%s.npy' % run_label), i1_scaled_iq)

    return result_df, min_scaled_iq, i1_scaled_iq, goal_iq


def interp_iq(data_iq, goal, q_max=0.2):

    q_grid = goal[:,0]
    data_int = np.zeros((len(q_grid), data_iq.shape[1]))
    data_int[:, 0] = q_grid

    # interpolate calculated data to be on the intended grid
    if data_iq.shape[1] == 3:
        print ('BE AWARE: replacing 3rd column in calculated I(Q) with '
               'goal error')
        interp_data = interpolate.splrep(data_iq[:, 0], data_iq[:, 1],
                                         w=1.0/data_iq[:,2])
        data_int[:, 1] = interpolate.splev(q_grid, interp_data)
        data_int[:, 2] = goal[:, 2]
    else:
        print ('Interpolating %d calculated I(Q) curves' %
               (data_iq.shape[1]-1))
        # todo: find a way to do this w/o looping
        for i in xrange(1, data_iq.shape[1]):
            interp_data = interpolate.splrep(data_iq[:, 0], data_iq[:, i])
            interp_iq = interpolate.splev(q_grid, interp_data)
            data_int[:, i] = interp_iq #/ interp_iq[0] # normalize to 1

    data_iq = data_int

    return data_iq


def polyspace(x1, x2, p, n):
    '''
    Usage:
       poly = polyspace(x1, x2, P, N)

    Purpose:
       create an array parabolicly spaced between x1 and x2

    Parameter(s):
       x1: first value of the array
       x2: second value of the array
       P:  polynomial power

    Return(s):
       poly: array of quadratically spaced values
    '''
    grid = np.linspace(x1 ** (1.0 / p), x2 ** (1.0 / p), n) ** p
    grid[0], grid[-1] = x1, x2
    return grid


def get_r(rf_data, mt_data):
    r, _ = get_r_components(rf_data, mt_data)

    return r


def get_r_components(rf_data, mt_data):
    # R value as defined by doi: 10.1042/bj2670203
    diff = np.abs(mt_data[:, 1] - rf_data[:, 1])
    norm = np.abs(rf_data[:, 1]).sum()
    components = diff / norm
    r = components.sum()

    return r, components


def get_f(rf_data, mt_data):
    f, _ = get_f2_components(rf_data, mt_data)

    return f

def get_f2_components(rf_data, mt_data):
    # F-factor as defined by doi: 10.1016/S0006-3495(98)77984-6
    diff = (np.log(mt_data[:, 1]) - np.log(rf_data[:, 1]))**2
    n = len(mt_data[:, 1])
    f2 = diff / n
    f = np.sqrt(np.nansum(f2)) # log produces nan from negative data

    return f, f2


def get_wr(rf_data, mt_data):
    wr, _ = get_wr_components(rf_data, mt_data)

    return wr


def get_wr_components(rf_data, mt_data):
    # weighted R-factor, correctly defined in doi: 10.1107/S0021889803000220
    w = rf_data[:, 0]
    w2 = w**2
    i1 = mt_data[:, 1]
    i2 = rf_data[:, 1]
    num = (w * (i1 - i2))**2
    den = w2 * i1 * i2
    wr_components = num / den.sum()
    wr = wr_components.sum()

    return wr, wr_components

def get_x2_norm(vals):
    '''
    calculate X^2 based on the assumption that vals are sampled from a normal
    distribution (taken from John R. Taylor "An Introduction to Error
    Analysis", pgs )
    '''
    # calculate the mean, standard deviation, and n
    # (each of these removes a dof)
    mean = np.mean(vals)
    std = np.std(vals)
    n = len(vals)

    # calculate the x2
    edges = np.linspace(mean - 2 * std, mean + 2 * std, 5)
    obs = np.arange(1, 9, dtype=np.float).reshape((2,4)).transpose()
    obs[:,1], _ = np.histogram(vals, bins=edges)
    exp = np.arange(1,13, dtype=np.float).reshape((3,4)).transpose()
    exp[:,1] = np.array([0.16, 0.34, 0.34, 0.16]) * n
    exp[:,2] = np.sqrt(exp[:,1])
    x2 = get_x2(exp, obs, dof=1)

    # get the probability of that distribution
    p = [[0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 6],
         [100, 62, 48, 39, 32, 26, 22 ,19, 16, 8, 5, 3, 1]]
    p = np.array(p).transpose()
    i_p_vals = np.argmin(np.abs(p[:,0]-x2))

    return x2, p[i_p_vals, 1], mean, std


def get_x2(rf_data, mt_data, dof=None):
    x2, _ = get_x2_components(rf_data, mt_data, dof=dof)

    return x2


def get_x2_components(rf_data, mt_data, dof=None):
    diff = mt_data[:, 1] - rf_data[:, 1]
    diff2 = diff * diff
    er2 = rf_data[:, 2] * rf_data[:, 2]
    if not dof:
        dof = len(rf_data)
    components = (diff2 / er2) / dof
    x2 = components.sum()

    return x2, components


def match_poly(in_data, rf_data):
    """
    determine the scale and offset to match the input data to the
    reference data using a polynomial fit

    Parameters
    ----------
    in_data:
        input data to match to the rf_data (should be Nx2 np.array)
    rf_data:
        reference data for matching the in_data (should be Nx3 np.array)

    Returns
    -------
    mt_data: version of in_data matched to the reference data
    scale:   scale factor applied to the input data
    offset:  offset applied to the input data

    Notes
    --------
    has the option to use error bars to weight the matching

    See also
    --------
    match_lstsq, scale, scale_offset

    """
    assert (in_data[:, 0] - rf_data[:, 0]).sum() == 0, ('mismatch between input and'
                                                        ' reference x-grid')
    try:
        weights = 1 / rf_data[:, 2]
    except:
        weights = None
    offset, scale = np.polynomial.polynomial.polyfit(
        in_data[:, 1], rf_data[:, 1], 1, w=weights)
    mt_data = np.vstack([in_data[:, 0], scale * in_data[:, 1] + offset]).T

    x2 = get_x2(rf_data, mt_data)

    return mt_data, scale, offset, x2


def scale_i0(in_data, rf_data):
    """
    determine the scale to match the I(0) or first value in the reference
    and input data

    Parameters
    ----------
    in_data:
        input data to match to the rf_data (should be Nx2 np.array)
    rf_data:
        reference data for matching the in_data (should be Nx3 np.array)

    Returns
    -------
    mt_data: version of in_data matched to the reference data
    scale:   scale factor applied to the input data
    x2:      X^2 comparison between the reference data and matched input data

    See also
    --------
    match_poly, match_lstsq, scale_offset, scale

    """
    assert (in_data[:, 0] - rf_data[:, 0]).sum() == 0, (
        'mismatch between input and reference x-grid')

    scale = rf_data[0, 1] / in_data[0, 1]

    mt_data = np.vstack([in_data[:, 0], scale * in_data[:, 1]]).T

    R = get_r(rf_data, mt_data)

    return mt_data, scale, R


def scale(in_data, rf_data):
    """
    determine the scale to match the input data to the reference
    data by minimizing the x2 calculation (critical that the
    error estimates are reasonable for all Q values)

    Parameters
    ----------
    in_data:
        input data to match to the rf_data (should be Nx2 np.array)
    rf_data:
        reference data for matching the in_data (should be Nx3 np.array)

    Returns
    -------
    mt_data: version of in_data matched to the reference data
    scale:   scale factor applied to the input data
    x2:      X^2 comparison between the reference data and matched input data

    See also
    --------
    match_poly, match_lstsq, scale_offset

    """
    assert (in_data[:, 0] - rf_data[:, 0]).sum() == 0, (
        'mismatch between input and reference x-grid')

    sigma2 = rf_data[:, 2] * rf_data[:, 2]
    scale = ((rf_data[:, 1] * in_data[:, 1] / sigma2).sum() /
             (in_data[:, 1] * in_data[:, 1] / sigma2).sum())

    mt_data = np.vstack([in_data[:, 0], scale * in_data[:, 1]]).T

    x2 = get_x2(rf_data, mt_data)

    return mt_data, scale, x2


def offset(in_data, rf_data):
    """
    determine the offset to match the input data to the reference
    data by minimizing the x2 calculation

    Parameters
    ----------
    in_data:
        input data to match to the rf_data (should be Nx2 np.array)
    rf_data:
        reference data for matching the in_data (should be Nx3 np.array)

    Returns
    -------
    mt_data: version of in_data matched to the reference data
    offset:   offset applied to the input data
    x2:      X^2 comparison between the reference data and matched input data

    See also
    --------
    match_poly, match_lstsq, scale_offset, scale

    """
    assert (in_data[:, 0] - rf_data[:, 0]).sum() == 0, ('mismatch between input and'
                                                        ' reference x-grid')

    sigma2 = rf_data[:, 2] * rf_data[:, 2]
    a = (rf_data[:, 1] / sigma2).sum()
    b = (in_data[:, 1] / sigma2).sum()
    c = (1 / sigma2).sum()
    offset = (a - b) / c

    mt_data = np.vstack([in_data[:, 0], in_data[:, 1] + offset]).T

    x2 = get_x2(rf_data, mt_data)

    return mt_data, offset, x2


def scale_offset(in_data, rf_data):
    """
    determine the scale and offset to match the input data to the reference
    data by minimizing the x2 calculation
    \chi^2 = \frac{1}{N_q-1}\sum_{i=1}^{N_q} \frac{\left[cI_{s_i}(Q) + I_c - I_{e_i}(Q)\right]^2}{\sigma_i^2(Q)}

    Parameters
    ----------
    in_data:
        input data to match to the rf_data (should be Nx2 np.array)
    rf_data:
        reference data for matching the in_data (should be Nx3 np.array)

    Returns
    -------
    mt_data: version of in_data matched to the reference data
    scale:   scale factor applied to the input data
    offset:  offset applied to the input data
    x2:      X^2 comparison between the reference data and matched input data

    See also
    --------
    match_poly, match_lstsq, scale

    """
    small = 1E-4  # small parameter
    assert np.allclose((in_data[:, 0] - rf_data[:, 0]).sum(), 0, atol=small), (
        'mismatch between input and reference x-grid')

    sigma2 = rf_data[:, 2] * rf_data[:, 2]
    a = (rf_data[:, 1] / sigma2).sum()
    b = (in_data[:, 1] / sigma2).sum()
    c = (1 / sigma2).sum()
    d = (rf_data[:, 1] * in_data[:, 1] / sigma2).sum()
    e = (in_data[:, 1] * in_data[:, 1] / sigma2).sum()

    offset = (a * e - b * d) / (c * e - b * b)
    scale = (c * d - b * a) / (c * e - b * b)

    mt_data = np.vstack([in_data[:, 0], scale * in_data[:, 1] + offset]).T

    x2 = get_x2(rf_data, mt_data)

    return mt_data, scale, offset, x2


def match_lstsq(in_data, rf_data):
    """
    determine the scale and offset to match the input data to the reference
    data using a lstsq fit

    Parameters
    ----------
    in_data:
        input data to match to the rf_data (should be Nx2 np.array)
    rf_data:
        reference data for matching the in_data (should be Nx2 np.array)

    Returns
    -------
    mt_data: version of in_data matched to the reference data
    scale:   scale factor applied to the input data
    offset:  offset applied to the input data

    Notes
    --------
    does not use error bars to weight the data

    See also
    --------
    match_poly, scale_offset, scale

    """
    assert (in_data[:, 0] - rf_data[:, 0]).sum() == 0, ('mismatch between input and'
                                                        ' reference x-grid')

    # could implement weights by changeing the second column to be 1/error
    A = np.vstack([in_data[:, 1], np.ones(len(in_data))]).T
    scale, offset = np.linalg.lstsq(A, rf_data[:, 1])[0]
    mt_data = np.vstack([in_data[:, 0], scale * in_data[:, 1] + offset]).T

    x2 = get_x2(rf_data, mt_data)

    return mt_data, scale, offset, x2


def scale_wr(in_data, rf_data):
    """
    determine the scale to use to match the input data to the reference
    data by minimizing the weighted R-factor

    Parameters
    ----------
    in_data:
        input data to match to the rf_data (should be Nx2 np.array)
    rf_data:
        reference data for matching the in_data (should be Nx2 np.array)

    Returns
    -------
    mt_data: version of in_data matched to the reference data
    scale:   scale factor applied to the input data
    offset:  offset applied to the input data

    Notes
    --------
    does not use error bars to weight the data

    See also
    --------
    match_poly, scale_offset, scale, match_lstsq

    """
    assert (in_data[:, 0] - rf_data[:, 0]).sum() == 0, ('mismatch between input and'
                                                        ' reference x-grid')

    w = rf_data[:, 0]
    w2 = w**2
    i1 = rf_data[:, 1]
    i2 = in_data[:, 1]
    # s1 = np.sum(i1 * i2 * w2) / np.sum(i2 * w2) # bogus units from publication
    scale = np.sqrt(np.sum(i1**2 * w2) / np.sum(i2**2 * w2)) # actual minimum

    mt_data = np.copy(in_data)
    mt_data[:,1:] *= scale
    wr = get_wr(rf_data, mt_data)

    return mt_data, scale, wr


def kratky(iq_data):
    '''
    iq_data should be a multi-column 2-D numpy array
    first column: Q
    returns a numpy array the same size where columns >= 2 are muliplied by q^2
    '''
    q = np.diag(iq_data[:, 0])
    q2 = q.dot(q)
    kratky_data = np.zeros(iq_data.shape)
    kratky_data[:, 0] = iq_data[:, 0]
    kratky_data[:, 1:] = q2.dot(iq_data[:, 1:])

    return kratky_data


def compare_match(rf_data, in_data, dummy_scale=None, dummy_offset=None,
                  log=False):
    scale = np.zeros(4)
    offset = np.zeros(4)
    x2 = np.zeros(4)
    mt_data_polyf, scale[0], offset[0], x2[0] = match_poly(in_data, rf_data)
    mt_data_lstsq, scale[1], offset[1], x2[1] = match_lstsq(in_data, rf_data)
    mt_data_scale, scale[2], x2[2] = scale(in_data, rf_data)
    mt_data_sclof, scale[3], offset[3], x2[3] = scale_offset(in_data, rf_data)
    info = []
    for i in xrange(4):
        info.append('s=%0.1f, o=%0.3f, x2=%0.1f' % (scale[i], offset[i], x2[i]))
    if dummy_scale and dummy_offset:
        ref_str = ', s=%0.1f, o=%0.1f' % (dummy_scale, dummy_offset)
    else:
        ref_str = ''

    fig = plt.figure()
    plt.subplot(2, 1, 1)
    plt.errorbar(rf_data[:, 0], rf_data[:, 1], yerr=rf_data[:, 2], fmt='o',
                 label='reference' + ref_str)
    plt.plot(mt_data_polyf[:, 0], mt_data_polyf[:, 1], 's',
             label='poly, %s' % info[0])
    plt.plot(mt_data_sclof[:, 0], mt_data_sclof[:, 1], '^',
             label='sc-of, %s' % info[3])
    plt.plot(mt_data_lstsq[:, 0], mt_data_lstsq[:, 1], '<',
             label='lstsq, %s' % info[1])
    plt.plot(mt_data_scale[:, 0], mt_data_scale[:, 1], '>',
             label='scale, %s' % info[2])
    plt.plot(in_data[:, 0], in_data[:, 1], 'o', label='input')
    if log:
        plt.yscale('log')
        plt.xscale('log')
        plt.axis('tight')
        plt.legend(loc=3)
    else:
        plt.legend(loc=1)

    tl1 = plt.title("Comparison of different match methods")

    plt.subplot(2, 1, 2)
    ref = rf_data[:, 1]
    plt.errorbar(rf_data[:, 0], rf_data[:, 1] - ref, yerr=rf_data[:, 2], fmt='o',
                 label='reference')
    plt.plot(mt_data_polyf[:, 0], mt_data_polyf[:, 1] - ref, 's', label='poly')
    plt.plot(mt_data_sclof[:, 0], mt_data_sclof[
             :, 1] - ref, '^', label='sc-of')
    plt.plot(mt_data_lstsq[:, 0], mt_data_lstsq[
             :, 1] - ref, '<', label='lstsq')
    plt.plot(mt_data_scale[:, 0], mt_data_scale[
             :, 1] - ref, '>', label='scale')
    tl2 = plt.title('residuals')
    if log:
        plt.xscale('log')
        plt.axis('tight')
    # plt.show()

    # i_sim = data_iq[1:,1]
    # i_exp = goal_iq[1:,1]
    # sigma = goal_iq[1:,2]
    # cs = np.linspace(0,20,201)
    # x2 = np.zeros(cs.shape)
    # dx2 = np.zeros(cs.shape)
    # for i, c in enumerate(cs):
    # x2[i] = ((c*i_sim-i_exp)**2/sigma**2).sum()/i_sim.shape[0]
    # dx2[i] = (2*I_sim*(c*i_sim-i_exp)/sigma**2).sum()/i_sim.shape[0]

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.hold()
    # plt.plot(cs, x2, label='x2')
    # plt.plot(cs, dx2, label='dx2/dc')
    # label_str = 'c=%0.3f' % c0
    # plt.axvline(x=c0, color='r', label=label_str)
    # plt.plot(cs,np.zeros(cs.shape),'k')
    # lg = plt.legend()
    # lg.draw_frame(False)
    # plt.title('best scale factor')
    # plt.show()


def examine_rg_i0(do_plot=False):
    df = load_rg_csv()

    di = ['diAtek010c050']
    tri = ['triEtek010c050', 'triEtek050c050',
           'triEtek100c050', 'triEtek010Mg1c050']
    tet = ['tetraAtek010c050', 'tetraAtek050c050',
           'tetraAtek100c050', 'tetraAtek010Mg1c050']
    data_files = {'di': di, 'tri': tri, 'tet': tet}

    if do_plot:
        plt.figure()
        x_range = [0, 120]
        #### Rg Subplots ####
        plt.subplot(3, 2, 1)
        plt.plot(df['KCl'].loc[di], df['Rg grp'].loc[di], 'o', label='GNOM rp',
                 markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[di], df['Rg grl'].loc[di],
                     df['RgEr grl'].loc[di], fmt='s', label='GNOM rl',
                     markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[di], df['Rg'].loc[di],
                     df['RgEr'].loc[di], fmt='>', label='MATLAB',
                     markeredgecolor='none')
        lg = plt.legend()
        plt.ylabel('Rg')
        plt.title('Dimer Rg comparison')
        plt.xlim(x_range)

        plt.subplot(3, 2, 3)
        plt.plot(df['KCl'].loc[tri], df['Rg grp'].loc[tri],
                 'o', label='GNOM rc', markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[tri], df['Rg grl'].loc[tri],
                     df['RgEr grl'].loc[tri], fmt='s', label='GNOM rl',
                     markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[tri], df['Rg'].loc[tri],
                     df['RgEr'].loc[tri], fmt='>', label='MATLAB',
                     markeredgecolor='none')
        # lg = plt.legend()
        plt.title('Trimer Rg comparison')
        plt.ylabel('Rg')
        plt.xlim(x_range)

        plt.subplot(3, 2, 5)
        plt.plot(df['KCl'].loc[tet], df['Rg grp'].loc[tet],
                 'o', label='GNOM rc', markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[tet], df['Rg grl'].loc[tet],
                     df['RgEr grl'].loc[tet], fmt='s', label='GNOM rl',
                     markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[tet], df['Rg'].loc[tet],
                     df['RgEr'].loc[tet], fmt='>', label='MATLAB',
                     markeredgecolor='none')
        # lg = plt.legend()
        plt.title('Tetramer Rg comparison')
        plt.ylabel('Rg')
        plt.xlabel(r'[KCl]')
        plt.xlim(x_range)

        #### I(0) Subplots ####
        plt.subplot(3, 2, 2)
        plt.plot(df['KCl'].loc[di], df['I0 grp'].loc[di],
                 'o', label='GNOM rc', markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[di], df['I0 grl'].loc[di],
                     df['I0Er grl'].loc[di], fmt='s', label='GNOM rl',
                     markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[di], df['I0'].loc[di],
                     df['I0Er'].loc[di], fmt='>', label='MATLAB',
                     markeredgecolor='none')
        plt.title('Dimer I(0) comparison')
        plt.ylabel('I(0)')
        plt.xlim(x_range)

        plt.subplot(3, 2, 4)
        plt.plot(df['KCl'].loc[tri], df['I0 grp'].loc[tri],
                 'o', label='GNOM rc', markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[tri], df['I0 grl'].loc[tri],
                     df['I0Er grl'].loc[tri], fmt='s', label='GNOM rl',
                     markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[tri], df['I0'].loc[tri],
                     df['I0Er'].loc[tri], fmt='>', label='MATLAB',
                     markeredgecolor='none')
        plt.title('Dimer I(0) comparison')
        plt.ylabel(r'I(0)')
        plt.xlim(x_range)

        plt.subplot(3, 2, 6)
        plt.plot(df['KCl'].loc[tet], df['I0 grp'].loc[tet],
                 'o', label='GNOM rc', markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[tet], df['I0 grl'].loc[tet],
                     df['I0Er grl'].loc[tet], fmt='s', label='GNOM rl',
                     markeredgecolor='none')
        plt.errorbar(df['KCl'].loc[tet], df['I0'].loc[tet],
                     df['I0Er'].loc[tet], fmt='>', label='MATLAB',
                     markeredgecolor='none')
        plt.title('Dimer I(0) comparison')
        plt.ylabel(r'I(0)')
        plt.xlabel(r'[KCl]')
        plt.xlim(x_range)

        # plt.show()

    return df, data_files


def fig_rg_v_conc():
    df = load_rg_csv()
    series = []
    labels = []

    series.append(['N12merH5Mg1b_zeroCon', 'N12merH5Mg1bx4', 'N12merH5Mg1bx2',
                   'N12merH5Mg1bx1'])
    labels.append(r'12x167 H5 1mM $Mg^{2+}$')

    series.append(['N12merH5TE_zeroCon', 'N12merH5TEx4', 'N12merH5TEx2',
                   'N12merH5TEx1'])
    labels.append(r'12x167 H5 10mM $K^+$')

    series.append(
        ['N12merH5Mg1a_zeroCon', 'N12merH5Mg1ax4', 'N12merH5Mg1ax2', 'N12merH5Mg1ax1'])
    labels.append(r'12x167 H5 1mM $Mg^{2+}$')

    series.append(
        ['N12merTE_zeroCon', 'N12merTEx4', 'N12merTEx2', 'N12merTEx1'])
    labels.append(r'12x167 10mM $K^{+}$')

    series.append(['N12merMg1_zeroCon', 'N12merMg1x4', 'N12merMg1x2',
                   'N12merMg1x1'])
    labels.append(r'12x167 1mM $Mg^{2+}$')

    series.append(['N4merH5TE_zeroCon', 'N4merH5TEx4', 'N4merH5TEx2',
                   'N4merH5TEx1'])
    labels.append(r'4x167 H5 10mM $K^{+}$')

    series.append(['N4merH5Mg1_zeroCon', 'N4merH5Mg1x4', 'N4merH5Mg1x2',
                   'N4merH5Mg1x1'])
    labels.append(r'4x167 H5 1mM $Mg^{2+}$')

    series.append(['N4merTE_zeroCon', 'N4merTEx4', 'N4merTEx2', 'N4merTEx1'])
    labels.append(r'4x167 10mM $K^{+}$ c')

    series.append(
        ['N4merMg1_zeroCon', 'N4merMg1x4', 'N4merMg1x2', 'N4merMg1x1'])
    labels.append(r'4x167 1mM $Mg^{2+}$ c')

    series.append(['GW3merTek010_zeroCon', 'GW3merAtek010x8', 'GW3merAtek010x4',
                   'GW3merAtek010x2', 'GW3merAtek010x1'])
    labels.append(r'3x167 10mM $K^+$ c')

    series.append(['GW3merTek050_zeroCon', 'GW3merAtek050x4', 'GW3merAtek050x2',
                   'GW3merAtek050x1'])
    labels.append(r'3x167 50mM $K^+$ c')

    series.append(['GW3merTek100_zeroCon', 'GW3merAtek100x4', 'GW3merAtek100x2',
                   'GW3merAtek100x1'])
    labels.append(r'3x167 100mM $K^+$ c')

    series.append(['GW3merTek200_zeroCon', 'GW3merAtek200x4', 'GW3merAtek200x2',
                   'GW3merAtek200x1'])
    labels.append(r'3x167 200mM $K^+$ c')

    series.append(['diAtek010_zeroCon', 'diAtek010c012', 'diAtek010c025',
                   'diAtek010c050'])
    labels.append(r'2x167 10mM $K^+$ na')

    series.append(['diBtek010_zeroCon', 'diBtek010c012', 'diBtek010c025',
                   'diBtek010c050'])
    labels.append(r'2x167 10mM $K^+$ nb')

    series.append(['triDtek010_zeroCon', 'triDtek010c012', 'triDtek010c025',
                   'triDtek010c050'])
    labels.append(r'3x167 10mM $K^+$ nd')

    series.append(['triEtek010_zeroCon', 'triEtek010c012', 'triEtek010c025',
                   'triEtek010c050'])
    labels.append(r'3x167 10mM $K^+$ ne')

    series.append(['triEtek050_zeroCon', 'triEtek050c012', 'triEtek050c025',
                   'triEtek050c050'])
    labels.append(r'3x167 50mM $K^+$ ne')

    series.append(['triEtek100_zeroCon', 'triEtek100c012', 'triEtek100c025',
                   'triEtek100c050'])
    labels.append(r'3x167 100mM $K^+$ ne')

    series.append(['triEMg1_zeroCon', 'triEtek010Mg1c012', 'triEtek010Mg1c025',
                   'triEtek010Mg1c050'])
    labels.append(r'3x167 1mM $Mg^{2+}$ ne')

    series.append(['triFtek010_zeroCon', 'triFtek010c012', 'triFtek010c025',
                   'triFtek010c050'])
    labels.append(r'3x167 10mM $K^+$ nf')

    series.append(['triFtek050_zeroCon', 'triFtek050c012', 'triFtek050c025',
                   'triFtek050c050'])
    labels.append(r'3x167 50mM $K^+$ nf')

    series.append(['tetraAtek010_zeroCon', 'tetraAtek010c012a', 'tetraAtek010c025',
                   'tetraAtek010c050'])
    labels.append(r'4x167 10mM $K^+$ na')

    series.append(['tetraAtek050_zeroCon', 'tetraAtek050c012', 'tetraAtek050c025',
                   'tetraAtek050c050'])
    labels.append(r'4x167 50mM $K^+$ na')

    series.append(['tetraAtek100_zeroCon', 'tetraAtek100c012', 'tetraAtek100c025',
                   'tetraAtek100c050'])
    labels.append(r'4x167 100mM $K^+$ na')

    series.append(['tetraAMg1_zeroCon', 'tetraAtek010Mg1c025',
                   'tetraAtek010Mg1c050', 'tetraAtek010Mg1c012'])
    labels.append(r'4x167 1mM $Mg^{2+}$ na')

    series.append(
        ['tetraCtek010_zeroCon', 'tetraCtek010c012', 'tetraCtek010c050'])
    labels.append(r'4x167 10mM $K^+$ nc')

    series.append(['tetraCtek050_zeroCon', 'tetraCtek050c012', 'tetraCtek050c025',
                   'tetraCtek050c050'])
    labels.append(r'4x167 50mM $K^+$ nc')

    plt.figure()
    x_range = [-1, 1]

    for i in xrange(len(series)):
        plt.errorbar(df['conc'].loc[series[i]], df['Rg'].loc[series[i]],
                     df['RgEr'].loc[series[i]], label=labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i, '--'),
                     mec=gp.qual_color(i), mfc='none', ms=15)

    lg = plt.legend(loc=0, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'mg/mL')
    plt.title(r'$R_g$ vs mg/mL comparison')
    plt.xlim(x_range)
    # plt.show()

    return


def fig_sub_i0_v_conc(show=False):
    df = load_rg_csv()
    chess = []
    nsls = []
    chess_labels = []
    nsls_labels = []
    l_2x167 = []
    l_3x167 = []
    l_4x167 = []
    l_gh5 = []

    # chess.append(['c000_12x167_k010',
                  # 'c125_12x167_k010',
                  # 'c250_12x167_k010',
                  # 'c500_12x167_k010'])
    # chess_labels.append(r'12x167: 10 mM $K^{+}$')

    # chess.append(['c000_12x167_mg1',
                  # 'c125_12x167_mg1',
                  # 'c250_12x167_mg1',
                  # 'c500_12x167_mg1'])
    # chess_labels.append(r'12x167: 1 mM $Mg^{2+}$')

    # chess.append(['c000_12x167_h5_k010',
                  # 'c125_12x167_h5_k010',
                  # 'c250_12x167_h5_k010',
                  # 'c500_12x167_h5_k010'])
    # chess_labels.append(r'12x167 gH5: 10 mM $K^+$')

    # chess.append(['c000_12x167_h5_mg1',
                  # 'c125_12x167_h5_mg1',
                  # 'c250_12x167_h5_mg1',
                  # 'c500_12x167_h5_mg1'])
    # chess_labels.append(r'12x167 gH5: 1 mM $Mg^{2+}$')

    tmp_list = ['c000_3x167_k010',
                  # 'c068_3x167_k010',
                  'c125_3x167_k010',
                  'c250_3x167_k010',
                  'c500_3x167_k010']
    chess.append(tmp_list)
    l_3x167 += tmp_list[1:]
    chess_labels.append(r'3x167: 10 mM $K^+$')

    tmp_list = ['c000_3x167_k050',
                'c125_3x167_k050',
                'c250_3x167_k050',
                'c500_3x167_k050']

    chess.append(tmp_list)
    l_3x167 += tmp_list[1:]
    chess_labels.append(r'3x167: 50 mM $K^+$')

    tmp_list = ['c000_3x167_k100',
                'c125_3x167_k100',
                'c250_3x167_k100',
                'c500_3x167_k100']
    chess.append(tmp_list)
    l_3x167 += tmp_list[1:]
    chess_labels.append(r'3x167: 100 mM $K^+$')

    tmp_list = ['c000_3x167_k200',
                'c125_3x167_k200',
                'c250_3x167_k200',
                'c500_3x167_k200']
    chess.append(tmp_list)
    l_3x167 += tmp_list[1:]
    chess_labels.append(r'3x167: 200 mM $K^+$')

    tmp_list = ['c000_4x167_k010',
                'c125_4x167_k010',
                'c250_4x167_k010',
                'c500_4x167_k010']
    nsls.append(tmp_list)
    l_4x167 += tmp_list[1:]
    nsls_labels.append(r'4x167: 10 mM $K^+$')

    tmp_list = ['c000_4x167_k050',
                'c125_4x167_k050',
                'c250_4x167_k050',
                'c500_4x167_k050']
    nsls.append(tmp_list)
    l_4x167 += tmp_list[1:]
    nsls_labels.append(r'4x167: 50 mM $K^+$')

    tmp_list = ['c000_4x167_k100',
                'c125_4x167_k100',
                'c250_4x167_k100',
                'c500_4x167_k100']
    nsls.append(tmp_list)
    l_4x167 += tmp_list[1:]
    nsls_labels.append(r'4x167: 100 mM $K^+$')

    tmp_list = ['c000_4x167_mg1',
                'c125_4x167_mg1',
                'c250_4x167_mg1',
                'c500_4x167_mg1']
    nsls.append(tmp_list)
    l_4x167 += tmp_list[1:]
    nsls_labels.append(r'4x167: 1 mM $Mg^{2+}$')

    chess.append(['c000_4x167_h5_k010',
                  'c200_4x167_h5_k010',
                  'c400_4x167_h5_k010',
                  'c800_4x167_h5_k010'])
    chess_labels.append(r'4x167 gH5: 10 mM $K^{+}$')

    chess.append(['c000_4x167_h5_mg1',
                  # 'c200_4x167_h5_mg1',
                  'c400_4x167_h5_mg1',
                  'c800_4x167_h5_mg1'])
    chess_labels.append(r'4x167 gH5: 1 mM $Mg^{2+}$')

    tmp_list = ['c000_2x167_k010',
                'c125_2x167_k010',
                'c250_2x167_k010',
                'c500_2x167_k010']
    nsls.append(tmp_list)
    l_2x167 += tmp_list[1:]
    nsls_labels.append(r'2x167: 10 mM $K^+$')

    df['I0_c'] = df['I0']/df['conc']
    df['I0Er_c'] = df['I0Er']/df['conc']

    x2_3x167, p_3x167, mean_3x167, std_3x167 = get_x2_norm(df.loc[l_3x167]['I0_c'])
    x2_4x167, p_4x167, mean_4x167, std_4x167 = get_x2_norm(df.loc[l_4x167]['I0_c'])
    print 'x2_3x167, p_3x167 =', x2_3x167, p_3x167
    print 'x2_4x167, p_4x167 =', x2_4x167, p_4x167

    fig = plt.figure(figsize=(6.335, 9))
    gs1 = GridSpec(2, 1)
    gs1.update(hspace=0)
    x_range = [0, 1.3]

    # SUBPLOT(1,1) a)
    ax = plt.subplot(gs1[0])
    for i in xrange(len(chess)):
        plt.errorbar(df['conc'].loc[chess[i]],
                     df['I0'].loc[chess[i]]/df['conc'].loc[chess[i]],
                     df['I0Er'].loc[chess[i]]/df['conc'].loc[chess[i]],
                     label=chess_labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i, '--'),
                     mec=gp.qual_color(i), mfc='none', ms=15, linewidth=2)
    plt.axhline(mean_3x167 + std_3x167, ls=':', c='DimGray')
    plt.axhline(mean_3x167 - std_3x167, ls=':', c='DimGray',
                label=r'$\sigma_{3\!\times 167}$')
    plt.ylabel(r'I(0)/(mg/mL)')
    plt.xlabel(r'mg/mL')
    # plt.ylim([0, 0.0009])
    plt.xlim(x_range)
    handles, labels = ax.get_legend_handles_labels()
    new_handles = handles[1:] + handles[:1]
    new_labels = labels[1:] + labels[:1]
    lg = ax.legend(new_handles, new_labels, loc='upper left', scatterpoints=1,
                   numpoints=1, bbox_to_anchor=(1, 1))
    lg.draw_frame(False)
    ax.get_xaxis().set_ticklabels([])
    ax.text(0.03, 0.92, r'(a) G1/CHESS', verticalalignment='bottom', fontweight='bold',
            horizontalalignment='left', transform=ax.transAxes)

    # SUBPLOT(2,1) b)
    ax = plt.subplot(gs1[1])
    for i in xrange(len(nsls)):
        plt.errorbar(df['conc'].loc[nsls[i]],
                     df['I0'].loc[nsls[i]]/df['conc'].loc[nsls[i]],
                     df['I0Er'].loc[nsls[i]]/df['conc'].loc[nsls[i]],
                     label=nsls_labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i, '--'),
                     mec=gp.qual_color(i), mfc='none', ms=15, linewidth=2)
    plt.axhline(mean_4x167 + std_4x167, ls=':', c='DimGray')
    plt.axhline(mean_4x167 - std_4x167, ls=':', c='DimGray',
                label=r'$\sigma_{4\!\times 167}$')
    handles, labels = ax.get_legend_handles_labels()
    new_handles = handles[1:] + handles[:1]
    new_labels = labels[1:] + labels[:1]
    lg = ax.legend(new_handles, new_labels, loc='upper left', scatterpoints=1,
                   numpoints=1, bbox_to_anchor=(1, 1))
    # lg = plt.legend(scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'I(0)/(mg/mL)')
    plt.xlabel(r'mg/mL')
    plt.ylim([140, 309])
    ax.set_yticks(ax.get_yticks()[:-1])
    plt.xlim(x_range)
    ax.text(0.03, 0.92, r'(b) X9/NSLS', verticalalignment='bottom', fontweight='bold',
            horizontalalignment='left', transform=ax.transAxes)

    fig.tight_layout()
    save_name = 'I0_v_mgmL'
    fig.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
    fig.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
    print 'view Rg vs mg/mL file: \nevince %s.eps &' % save_name
    if show:
        plt.show()
    else:
        plt.close()

    return



def fig_sub_rg_v_conc(show=False):
    df = load_rg_csv()
    dod = []
    tet = []
    tri = []
    di = []
    dod_labels = []
    tet_labels = []
    tri_labels = []
    di_labels = []

    dod.append(['c000_12x167_k010',
                'c125_12x167_k010',
                'c250_12x167_k010',
                'c500_12x167_k010'])
    dod_labels.append(r'10mM K$^{+}$')

    dod.append(['c000_12x167_mg1',
                'c125_12x167_mg1',
                'c250_12x167_mg1',
                'c500_12x167_mg1'])
    dod_labels.append(r'1mM Mg$^{2+}$')

    dod.append(['c000_12x167_h5_k010',
                'c125_12x167_h5_k010',
                'c250_12x167_h5_k010',
                'c500_12x167_h5_k010'])
    dod_labels.append(r'gH5 10mM K$^+$')

    dod.append(['c000_12x167_h5_mg1',
                'c125_12x167_h5_mg1',
                'c250_12x167_h5_mg1',
                'c500_12x167_h5_mg1'])
    dod_labels.append(r'gH5 1mM Mg$^{2+}$')

    tri.append(['c000_3x167_k010',
                # 'c068_3x167_k010',
                'c125_3x167_k010',
                'c250_3x167_k010',
                'c500_3x167_k010'])
    tri_labels.append(r'10mM K$^+$')

    tri.append(['c000_3x167_k050',
                'c125_3x167_k050',
                'c250_3x167_k050',
                'c500_3x167_k050'])
    tri_labels.append(r'50mM K$^+$')

    tri.append(['c000_3x167_k100',
                'c125_3x167_k100',
                'c250_3x167_k100',
                'c500_3x167_k100'])
    tri_labels.append(r'100mM K$^+$')

    tri.append(['c000_3x167_k200',
                'c125_3x167_k200',
                'c250_3x167_k200',
                'c500_3x167_k200'])
    tri_labels.append(r'200mM K$^+$')

    di.append(['c000_2x167_k010',
               'c125_2x167_k010',
               'c250_2x167_k010',
               'c500_2x167_k010'])
    di_labels.append(r'10mM K$^+$')

    tet.append(['c000_4x167_k010',
                'c125_4x167_k010',
                'c250_4x167_k010',
                'c500_4x167_k010'])
    tet_labels.append(r'10mM K$^+$')

    tet.append(['c000_4x167_k050',
                'c125_4x167_k050',
                'c250_4x167_k050',
                'c500_4x167_k050'])
    tet_labels.append(r'50mM K$^+$')

    tet.append(['c000_4x167_k100',
                'c125_4x167_k100',
                'c250_4x167_k100',
                'c500_4x167_k100'])
    tet_labels.append(r'100mM K$^+$')

    tet.append(['c000_4x167_mg1',
                'c125_4x167_mg1',
                'c250_4x167_mg1',
                'c500_4x167_mg1'])
    tet_labels.append(r'1mM Mg$^{2+}$')

    tet.append(['c000_4x167_h5_k010',
                'c200_4x167_h5_k010',
                'c400_4x167_h5_k010',
                'c800_4x167_h5_k010'])
    tet_labels.append(r'gH5 10mM K$^{+}$')

    tet.append(['c000_4x167_h5_mg1',
                # 'c200_4x167_h5_mg1',
                'c400_4x167_h5_mg1',
                'c800_4x167_h5_mg1'])
    tet_labels.append(r'gH5 1mM Mg$^{2+}$')

    fig = plt.figure(figsize=(12, 9))
    gs1 = GridSpec(2, 2)
    gs1.update(hspace=0)
    x_range = [-0.05, 1.5]

    # SUBPLOT(1,1) a)
    ax = plt.subplot(gs1[0])
    for i in xrange(len(tri)):
        plt.errorbar(df['conc'].loc[tri[i]], df['Rg'].loc[tri[i]],
                     df['RgEr'].loc[tri[i]], label=tri_labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i, '--'),
                     mec=gp.qual_color(i), mfc='none', ms=15, linewidth=2)

    plt.ylabel(r'$R_g$')
    plt.xlabel(r'mg/mL')
    # plt.title(r'3x167', x=0.3, y=0.92)
    # ylim = np.array(plt.ylim())
    # ylim[1] *= 1.08
    # plt.ylim(ylim)
    # tri_ylim = ylim
    tri_ylim = [90, 116]  # store for dimer ylim
    plt.ylim(tri_ylim)
    plt.xlim(x_range)
    lg = plt.legend(loc='center right', scatterpoints=1, numpoints=1)
    # lg = plt.legend(scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    ax.get_xaxis().set_ticklabels([])
    ax.text(0.03, 0.92, r'(a) 3x167', verticalalignment='bottom', fontweight='bold',
            horizontalalignment='left', transform=ax.transAxes)

    # SUBPLOT(2,1) b)
    ax = plt.subplot(gs1[2])
    for i in xrange(len(tet)):
        plt.errorbar(df['conc'].loc[tet[i]], df['Rg'].loc[tet[i]],
                     df['RgEr'].loc[tet[i]], label=tet_labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i, '--'),
                     mec=gp.qual_color(i), mfc='none', ms=15, linewidth=2)

    # lg = plt.legend(loc='upper right', scatterpoints=1, numpoints=1)
    lg = plt.legend(scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'mg/mL')
    # plt.title(r'4x167', x=0.3, y=0.92)
    # ylim = np.array(plt.ylim())
    # ylim[1] *= 1.03
    # plt.ylim(ylim)
    ax.text(0.03, 0.92, r'(b) 4x167', verticalalignment='bottom', fontweight='bold',
            horizontalalignment='left', transform=ax.transAxes)
    plt.xlim(x_range)
    plt.ylim([75, 150])
    ax.set_yticks(ax.get_yticks()[1:-1])


    x_range = [-0.05, 0.8]

    # SUBPLOT(1,2) c)
    ax = plt.subplot(gs1[1])
    # plt.figure()
    for i in xrange(len(di)):
        this_color = gp.qual_color(i)
        plt.errorbar(df['conc'].loc[di[i]], df['Rg'].loc[di[i]],
                     df['RgEr'].loc[di[i]], label=di_labels[i],
                     c=this_color, fmt=gp.symbol_order(i, '--'),
                     mec=this_color, mfc='none', ms=15, linewidth=2)

    lg = plt.legend(loc='center right', scatterpoints=1, numpoints=1)
    # lg = plt.legend(scatterpoints=1, numpoints=1)
    try:
        lg.draw_frame(False)
    except:
        pass
    plt.ylabel(r'$R_g$')
    # plt.xlabel(r'mg/mL')
    # plt.title('2x167', x=0.3, y=0.92)
    # plt.text(.2, .9, r'2x167: $R_g$ vs mg/mL comparison',
    # horizontalalignment='center',nn
    # transform=plt.transAxes)
    plt.xlim(x_range)
    plt.ylim(tri_ylim)
    ax.get_xaxis().set_ticklabels([])
    ax.text(0.03, 0.92, r'(c) 2x167', verticalalignment='bottom', fontweight='bold',
            horizontalalignment='left', transform=ax.transAxes)

    # SUBPLOT(2,2) d)
    ax = plt.subplot(gs1[3])
    for i in xrange(len(dod)):
        plt.errorbar(df['conc'].loc[dod[i]], df['Rg'].loc[dod[i]],
                     df['RgEr'].loc[dod[i]], label=dod_labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i, '--'),
                     mec=gp.qual_color(i), mfc='none', ms=15, linewidth=2)

    # lg = plt.legend(loc='upper right', scatterpoints=1, numpoints=1)
    lg = plt.legend(scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'mg/mL')
    # plt.title(r'12x167', x=0.3, y=0.92)
    ylim = np.array(plt.ylim())
    # ax.set_yticks(ax.get_yticks()[:-1])
    ylim[1] *= 1.03
    plt.ylim(ylim)
    plt.xlim(x_range)
    ax.text(0.03, 0.92, r'(d) 12x167', verticalalignment='bottom', fontweight='bold',
            horizontalalignment='left', transform=ax.transAxes)

    fig.tight_layout()
    save_name = 'Rg_v_mgmL'
    fig.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
    fig.savefig(save_name + '.eps', dpi=400, bbox_inches='tight')
    print 'view Rg vs mg/mL file: \nevince %s.eps &' % save_name
    if show:
        plt.show()
    else:
        plt.close()

    return


def fig_rg_v_salt(show=False):
    df = load_rg_csv()
    tri0 = ['c000_3x167_k010',
            'c000_3x167_k050',
            'c000_3x167_k100',
            'c000_3x167_k200']
    tri5 = ['c500_3x167_k010',
            'c500_3x167_k050',
            'c500_3x167_k100',
            'c500_3x167_k200']
    tet0 = ['c000_4x167_k010',
            'c000_4x167_k050',
            'c000_4x167_k100']
    tet5 = ['c500_4x167_k010',
            'c500_4x167_k050',
            'c500_4x167_k100']
    # triE = ['triEtek010_zeroCon', 'triEtek050_zeroCon',
    # 'triEtek100_zeroCon']
    di0 = ['c000_2x167_k010']
    di5 = ['c500_2x167_k010']
    # series = [di0, di5, tri0, tri5, tet0, tet5]
    # labels = ['2x167 (0.0 mg/mL)',
              # '2x167 (0.5 mg/mL)',
              # '3x167 (0.0 mg/mL)',
              # '3x167 (0.5 mg/mL)',
              # '4x167 (0.0 mg/mL)',
              # '4x167 (0.5 mg/mL)']
    series = [di0, tri0, tet0]
    labels = ['2x167 (0.0 mg/mL)',
              '3x167 (0.0 mg/mL)',
              '4x167 (0.0 mg/mL)']
    fig = plt.figure()
    x_range = [-10, 220]
    for i in xrange(len(series)):
        plt.errorbar(df['KCl'].loc[series[i]], df['Rg'].loc[series[i]],
                     df['RgEr'].loc[series[i]], label=labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i, '--'),
                     mec=gp.qual_color(i), mfc='none', ms=15)

    lg = plt.legend(loc=0, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'[$K^+$]')
    plt.title(r'$R_g$ vs [$K^+$] comparison')
    plt.xlim(x_range)
    if show:
        plt.show()
    fig.tight_layout()
    fig.savefig('Rg_v_salt.png')
    fig.savefig('Rg_v_salt.eps')
    plt.close()


def load_rg_csv():
    output_csv = ('/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/1406data'
                  '/rg_i0.csv')
    if op.exists(output_csv) and not debug:
        print '# READING Rg AND I(0) FROM CSV FILE'
        df = pd.read_csv(output_csv, sep='\t')
        # df = pd.read_csv(output_csv, sep=',')
        df.index = df['prefix']
    else:
        print '# READING Rg AND I(0) FROM MATLAB AND GNOM OUTPUT'
        df = combine_rg_i0()
        df.to_csv(output_csv, sep='\t')
    return df


def combine_rg_i0():
    data_dir = 'iqdata'

    # chess_dir = '/home/schowell/Dropbox/gw_phd/experiments/1406CHESS/'
    chess_dir = '/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/1406data/chess/'
    chess_matlab_file = 'rg_i0_chess.csv'
    chess_matlab_df = pd.read_csv(op.join(chess_dir, chess_matlab_file))
    chess_matlab_df = chess_matlab_df.set_index(
        chess_matlab_df['prefix']).drop('prefix', 1)

    # chess_gnom_files = glob.glob(op.join(chess_dir, data_dir, '*.out'))
    # chess_gnom_df = gnom_rg(chess_gnom_files)

    # nsls_dir = '/home/schowell/Dropbox/gw_phd/experiments/1406NSLS/'
    nsls_dir = '/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/1406data/nsls/'
    nsls_matlab_file = 'rg_i0_nsls.csv'
    nsls_matlab_df = pd.read_csv(op.join(nsls_dir, nsls_matlab_file))
    nsls_matlab_df = nsls_matlab_df.set_index(
        nsls_matlab_df['prefix']).drop('prefix', 1)

    # nsls_gnom_files = glob.glob(op.join(nsls_dir, data_dir, '*.out'))
    # nsls_gnom_df = gnom_rg(nsls_gnom_files)

    zero_dir = '/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/1406data/'
    zero_gnom_files = glob.glob(op.join(zero_dir, data_dir, '*.out'))
    zero_gnom_df = gnom_rg(zero_gnom_files)

    # chess_df = pd.concat([chess_gnom_df, chess_matlab_df])
    # nsls_df = pd.concat([nsls_gnom_df, nsls_matlab_df])
    # df = pd.concat([chess_df, nsls_df])
    combined_df = pd.concat([chess_matlab_df, nsls_matlab_df])
    combined_df.drop_duplicates(inplace=True)
    df = pd.concat([zero_gnom_df, combined_df], axis=1)
    df.index.name = 'prefix'

    KCl = []
    for prefix in df.index:
        if 'mg1' in prefix.lower():
            salt = 211
        elif 'k010' in prefix.lower() or 'TE' in prefix:
            salt = 10
        elif 'k100' in prefix.lower():
            salt = 100
        elif 'k200' in prefix.lower():
            salt = 200
        elif 'k050' in prefix.lower():
            salt = 50
        else:
            salt = np.nan
        if False:
            print '%s: ' % prefix, salt
        KCl.append(salt)
    df['KCl'] = KCl

    zero_df = df.loc[zero_gnom_df.index]
    zero_df['Ipr/Ig'] = zero_df['I0 gr']/zero_df['I0']
    zero_df['Rg_pr/Rg'] = zero_df['Rg gr']/zero_df['Rg']
    zero_df['dRg'] = np.abs(zero_df['Rg']-zero_df['Rg gr'])
    zero_df['dRg er'] = np.sqrt(zero_df['RgEr gr']**2 + zero_df['RgEr']**2)

    zero_output_csv = ('/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/'
                       '1406data/c000_rg_i0.csv')
    zero_df.to_csv(zero_output_csv, sep='\t')

    return df


def compare_iq(array_types, data_files, data_dir, data_ext, run_dirs,
                prefix='', do_plot=True, cutoff=None, best_dcd=False,
                fresh=False, run_label='', metric='r', N=100):
    all_result_dfs = []
    all_i1_data_iqs = []
    all_min_data_iqs = []
    all_goal_iqs = []
    all_data_files = []

    for array_type in array_types:
        for data_file in data_files[array_type]:
            rf_file = op.join(data_dir, data_file + data_ext) # reference file
            assert op.exists(rf_file), (
                "No such file: '%s'" % rf_file)

            df_list = []
            i1_data_iq_list = []
            min_data_iq_list = []
            for run_dir in run_dirs[array_type]:
                filter_dir = op.join(run_dir, data_file + '_filter')
                out_file = op.join(filter_dir, 'rg_x2%s.out' % run_label)
                if op.exists(out_file) and not fresh:
                    print 'loading rg and discrepancy from %s' % out_file
                    result_df = pd.read_csv(out_file, sep='\t')
                    i1_data_iq = np.load(op.join(filter_dir, 'i1_data_%s.npy' %
                                              run_label))
                    min_data_iq = np.load(op.join(filter_dir, 'min_data_%s.npy'
                                                  % run_label))
                    goal_iq = np.loadtxt(op.join(filter_dir, 'goal_%s.iq' %
                                                 run_label))
                else:
                    print ('loading rg then calculating discrepancy for %s' %
                           run_dir)
                    rf_data = np.loadtxt(rf_file)
                    rf_data[:, 1:] /= rf_data[0, 1] # normalize data
                    result_df, min_data_iq, i1_data_iq, goal_iq = compare_run_to_iq(
                        run_dir, rf_data, filter_dir, out_file=out_file,
                        metric=metric, run_label=run_label)
                run_name = run_dir.split('/')[-3] + '/' + run_dir.split('/')[-2]
                result_df['run'] = run_name
                df_list.append(result_df)
                i1_data_iq_list.append(i1_data_iq)
                min_data_iq_list.append(min_data_iq)

            # combine result DataFrames and data_iq arrays
            result_df = pd.concat(df_list)
            result_df.index = range(len(result_df))


            i1_key = 'i1_%s' % metric
            result_df.sort(i1_key, inplace=True)
            best200 = result_df.iloc[:200]
            best200.to_csv(data_file + '_%s_best.csv' % i1_key,
                           float_format='%5.10f', sep='\t')

            if metric == 'x2' or metric == 'wR':
                min_key = 'min_%s' % metric
                result_df.sort(min_key, inplace=True)
                best200 = result_df.iloc[:200]
                best200.to_csv(data_file + '_%s_best.csv' % min_key,
                           float_format='%5.10f', sep='\t')
            else:
                min_key = 'min_i0'

            q = i1_data_iq_list[0][:, :1]
            for (i, i1_data_iq) in enumerate(i1_data_iq_list):
                if i == 0:
                    i1_data_iqs = np.concatenate((q, i1_data_iq[:, 1:]), axis=1)
                else:
                    i1_data_iqs = np.concatenate((i1_data_iqs, i1_data_iq[:, 1:]),
                                                  axis=1)
                assert np.abs(np.sum(i1_data_iq[0, 1:] -goal_iq[0,1])) < 1e-5, (
                    'ERROR: Failed to scale data to I(0)')

            prefix = i1_key
            if run_label:
                prefix += '_%s' % run_label
            plot_discrepancy(result_df, i1_data_iqs, goal_iq, data_file,
                             prefix=prefix, key=i1_key,
                             residual=r'$R$-factor', N=N)
            plot_discrepancy(result_df, i1_data_iqs, goal_iq, data_file,
                             prefix=prefix, key=i1_key,
                             residual=r'$R$-factor', sub_label=False, N=N)

            for (i, min_data_iq) in enumerate(min_data_iq_list):
                if i == 0:
                    min_data_iqs = np.concatenate((q, min_data_iq[:, 1:]),
                                                  axis=1)
                else:
                    min_data_iqs = np.concatenate((min_data_iqs,
                                                   min_data_iq[:, 1:]),
                                                  axis=1)

            prefix = min_key
            if run_label:
                prefix += '_%s' % run_label
            plot_discrepancy(result_df, min_data_iqs, goal_iq, data_file,
                             prefix=prefix, key=min_key,
                             residual=r'$R$-factor', N=N)

            if cutoff or best_dcd:
                write_filter_output(run_dirs[array_type], df_list, cutoff,
                                    result_df, i1_data_iq_list, key=i1_key,
                                    best_dcd=best_dcd, data_file=data_file,
                                    label='_%s%s' % ( data_file, run_label),
                                    goal_iq=goal_iq, N=N)

            all_result_dfs.append(result_df)
            all_i1_data_iqs.append(i1_data_iqs)
            all_min_data_iqs.append(min_data_iqs)
            all_goal_iqs.append(goal_iq)
            all_data_files.append(data_file)

    return (all_result_dfs, all_i1_data_iqs, all_min_data_iqs, all_goal_iqs,
            all_data_files)



def compare_iq_all_metrics(array_types, data_files, data_dir, data_ext,
                           run_dirs, prefix='', do_plot=True, cutoff=None,
                           best_dcd=False, fresh=False, s=False, o=False):
    '''
    created: 14 Oct 2015
    purpose: compare the different discrepancy methods
    '''
    all_x2rg_dfs = []
    all_data_iqs = []
    all_goal_iqs = []
    all_data_files = []

    for array_type in array_types:
        for data_file in data_files[array_type]:
            rf_file = op.join(data_dir, data_file + data_ext)
            assert op.exists(rf_file), (
                "No such file: '%s'" % rf_file)

            df_list = []
            x2_data_iq_list = []
            wr_data_iq_list = []
            i1_data_iq_list = []
            for run_dir in run_dirs[array_type]:
                filter_dir = op.join(run_dir, data_file + '_filter')
                # run_label = '_s%do%d' % (s, o)
                run_label = ''
                out_file = op.join(filter_dir, 'rg_x2%s.out' % run_label)
                if op.exists(out_file) and not fresh:
                    print 'loading rg and discrepancy from %s' % out_file
                    result_df = pd.read_csv(out_file, sep='\t')
                    x2_data_iq = np.load(op.join(filter_dir, 'x2_data_iq%s.npy'
                                                 % run_label))
                    wr_data_iq = np.load(op.join(filter_dir, 'wr_data_iq%s.npy'
                                                 % run_label))
                    i1_data_iq = np.load(op.join(filter_dir, 'i1_data_iq%s.npy'
                                                 % run_label))
                    goal_iq = np.loadtxt(op.join(filter_dir, 'goal%s.iq' %
                                                 run_label))
                else:
                    print ('loading rg then calculating discrepancy for %s' %
                           run_dir)
                    rf_data = np.loadtxt(rf_file)
                    # rf_data[:, 1:] /= rf_data[0, 1] # normalize reference data
                    (result_df, x2_data_iq, wr_data_iq, i1_data_iq, goal_iq
                     )= compare_run_to_iq_all_metrics(run_dir, rf_data, filter_dir,
                                          out_file, run_label=run_label,
                                          match_range=[0.004, 0.2])
                run_name = run_dir.split('/')[-3] + '/' + run_dir.split('/')[-2]
                result_df['run'] = run_name
                df_list.append(result_df)
                x2_data_iq_list.append(x2_data_iq)
                wr_data_iq_list.append(wr_data_iq)
                i1_data_iq_list.append(i1_data_iq)
            if cutoff:
                write_filter_output(run_dirs[array_type], df_list, cutoff,
                                    best_dcd=best_dcd, label='_%s%s' %
                                    (data_file, run_label))

            # combine result DataFrames and data_iq arrays
            x2rg_df = pd.concat(df_list)
            x2rg_df.index = range(len(x2rg_df))

            # plot the residual values vs index when scaled using i1
            fig = plt.figure()
            gs1 = GridSpec(4, 4, hspace=0, wspace=0)#, left=0.1, right=0.9, bottom=0.075, top=0.925,

            ax0 = plt.subplot(gs1[0,:])
            ax0.plot(x2rg_df.i1_x2/x2rg_df.i1_x2.iloc[0], label='i1_x2')
            ax0.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax0.get_xaxis().set_ticks([])
            ax0.get_yaxis().set_ticks([])

            ax1 = plt.subplot(gs1[1,:])
            ax1.plot(x2rg_df.i1_wr/x2rg_df.i1_wr.iloc[0], label='i1_wr')
            ax1.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax1.get_xaxis().set_ticks([])
            ax1.get_yaxis().set_ticks([])

            ax2 = plt.subplot(gs1[2,:])
            ax2.plot(x2rg_df.i1_r/x2rg_df.i1_r.iloc[0], label='i1_r')
            ax2.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax2.get_xaxis().set_ticks([])
            ax2.get_yaxis().set_ticks([])

            ax3 = plt.subplot(gs1[3,:])
            ax3.plot(x2rg_df.i1_f/x2rg_df.i1_f.iloc[0], label='i1_f')
            ax3.legend(loc='upper left', bbox_to_anchor=(1, 1))
            plt.xlabel('Structure Number')
            ax3.get_xaxis().set_ticks([])

            save_name = data_file + '_i1_disc'
            fig.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
            plt.close()
            # plt.show()

            # plot the residual values vs index when scaled using high Q

            fig = plt.figure()
            gs2 = GridSpec(8, 4, hspace=0, wspace=0)

            ax0 = plt.subplot(gs2[0,:])
            ax0.plot(x2rg_df.mx2_x2/x2rg_df.mx2_x2.iloc[0], label='mx2_x2')
            ax0.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax0.get_xaxis().set_ticks([])
            ax0.get_yaxis().set_ticks([])

            ax1 = plt.subplot(gs2[1,:])
            ax1.plot(x2rg_df.mx2_wr/x2rg_df.mx2_wr.iloc[0], label='mx2_wr')
            ax1.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax1.get_xaxis().set_ticks([])
            ax1.get_yaxis().set_ticks([])

            ax2 = plt.subplot(gs2[2,:])
            ax2.plot(x2rg_df.mx2_r / x2rg_df.mx2_r.iloc[0],  label='mx2_r')
            ax2.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax2.get_xaxis().set_ticks([])
            ax2.get_yaxis().set_ticks([])

            ax3 = plt.subplot(gs2[3,:])
            ax3.plot(x2rg_df.mx2_f / x2rg_df.mx2_f.iloc[0],  label='mx2_f')
            ax3.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax3.get_xaxis().set_ticks([])
            ax3.get_yaxis().set_ticks([])

            ax4 = plt.subplot(gs2[4,:])
            ax4.plot(x2rg_df.mwr_x2/x2rg_df.mwr_x2.iloc[0], label='mwr_x2')
            ax4.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax4.get_xaxis().set_ticks([])
            ax4.get_yaxis().set_ticks([])

            ax5 = plt.subplot(gs2[5,:])
            ax5.plot(x2rg_df.mwr_wr/x2rg_df.mwr_wr.iloc[0], label='mwr_wr')
            ax5.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax5.get_xaxis().set_ticks([])
            ax5.get_yaxis().set_ticks([])

            ax6 = plt.subplot(gs2[6,:])
            ax6.plot(x2rg_df.mwr_r / x2rg_df.mwr_r.iloc[0],  label='mwr_r')
            ax6.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax6.get_xaxis().set_ticks([])
            ax6.get_yaxis().set_ticks([])

            ax7 = plt.subplot(gs2[7,:])
            ax7.plot(x2rg_df.mwr_f / x2rg_df.mwr_f.iloc[0],  label='mwr_f')
            ax7.legend(loc='upper left', bbox_to_anchor=(1, 1))
            ax7.get_xaxis().set_ticks([])
            ax7.get_yaxis().set_ticks([])

            plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
            plt.xlabel('Structure Number')

            save_name = data_file + '_hiQ_disc'
            fig.savefig(save_name + '.png', dpi=400, bbox_inches='tight')
            plt.close()
            # plt.show()


            x2rg_df.sort('i1_x2', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_i1_x2_best.csv',
                           float_format='%5.10f', sep='\t')

            x2rg_df.sort('i1_wr', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_i1_wr_best.csv',
                           float_format='%5.10f', sep='\t')

            x2rg_df.sort('i1_r', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_i1_r_best.csv',
                           float_format='%5.10f', sep='\t')

            x2rg_df.sort('i1_f', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_i1_f_best.csv',
                           float_format='%5.10f', sep='\t')

            x2rg_df.sort('mx2_x2', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_x2_best.csv',
                           float_format='%5.10f', sep='\t')

            x2rg_df.sort('mwr_wr', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_wr_best.csv',
                           float_format='%5.10f', sep='\t')

            x2rg_df.sort('mwr_r', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_mwr_r_best.csv',
                           float_format='%5.10f', sep='\t')

            x2rg_df.sort('mx2_r', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_mx2_r_best.csv',
                           float_format='%5.10f', sep='\t')


            x2rg_df.sort('mwr_f', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_mwr_f_best.csv',
                           float_format='%5.10f', sep='\t')

            x2rg_df.sort('mx2_f', inplace=True)
            best200 = x2rg_df.iloc[:200]
            best200.to_csv(data_file + '_mx2_f_best.csv',
                           float_format='%5.10f', sep='\t')

            q = x2_data_iq_list[0][:, :1]
            for (i, x2_data_iq) in enumerate(x2_data_iq_list):
                if i == 0:
                    x2_data_iqs = np.concatenate((q, x2_data_iq[:, 1:]), axis=1)
                else:
                    x2_data_iqs = np.concatenate((x2_data_iqs, x2_data_iq[:, 1:]),
                                                 axis=1)

            for (i, wr_data_iq) in enumerate(wr_data_iq_list):
                if i == 0:
                    wr_data_iqs = np.concatenate((q, wr_data_iq[:, 1:]), axis=1)
                else:
                    wr_data_iqs = np.concatenate((wr_data_iqs, wr_data_iq[:, 1:]),
                                                 axis=1)

            for (i, i1_data_iq) in enumerate(i1_data_iq_list):
                if i == 0:
                    i1_data_iqs = np.concatenate((q, i1_data_iq[:, 1:]), axis=1)
                else:
                    i1_data_iqs = np.concatenate((i1_data_iqs, i1_data_iq[:, 1:]),
                                                 axis=1)
                assert np.abs(np.sum(i1_data_iq[0, 1:] -goal_iq[0,1])) < 1e-5, (
                    'ERROR: Failed to scale data to I(0)')

            # i1 scale, x2
            prefix = 'i1_x2'
            plot_discrepancy(x2rg_df, i1_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'$\chi^2$', key='i1_x2')
            # x2 scale, x2
            prefix = 'mx2_x2'
            plot_discrepancy(x2rg_df, x2_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'$\chi^2$', key='mx2_x2')
            # i1 scale, wr
            prefix = 'i1_wr'
            plot_discrepancy(x2rg_df, i1_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'weighted $R$', key='i1_wr')
            # wr scale, wr
            prefix = 'mwr_wr'
            plot_discrepancy(x2rg_df, wr_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'weighted $R$', key='mwr_wr')
            # i1 scale, r
            prefix = 'i1_r'
            plot_discrepancy(x2rg_df, i1_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'$R$-factor', key='i1_r')
            # x2 scale, r
            prefix = 'mx2_r'
            plot_discrepancy(x2rg_df, x2_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'$R$-Factor', key='mx2_r')
            # wr scale, r
            prefix = 'mwr_r'
            plot_discrepancy(x2rg_df, wr_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'$R$-Factor', key='mwr_r')
            # i1 scale, f
            prefix = 'i1_f'
            plot_discrepancy(x2rg_df, i1_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'$F$-Factor', key='i1_f')
            # x2 scale, f
            prefix = 'mx2_f'
            plot_discrepancy(x2rg_df, x2_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'$F$-Factor', key='mx2_f')
            # wr scale, f
            prefix = 'mwr_f'
            plot_discrepancy(x2rg_df, wr_data_iqs, goal_iq, data_file,
                          prefix + run_label[1:],
                          residual=r'$F$-Factor', key='mwr_f')


def write_filter_output(run_dirs, df_list, cutoff, result_df, data_iq_list,
                        best_dcd=False, key='i1_wR', label='', data_file='',
                        catdcd_exe='/home/schowell/data/myPrograms/bin/catdcd',
                        goal_iq=[], do_align=False, N=100):
    txt_name = 'rglowweights%s.txt' % label
    index = 0

    with open(op.join(txt_name), 'w') as txt_file:
        txt_file.write('# file generated on FILL THIS IN\n')
        txt_file.write('# structure, Rg, weight\n')

        if best_dcd:
            # create a dcd with the best structures and extract best I(Q)s
            n_accept = 0
            best_iqs = []
            best_iqs.append(data_iq_list[0][:, 0]) # the Q-values

            # read a random pdb from a preceding directory, LIKELY TO FAIL
            random_pdb = glob.glob(op.split(op.split(
                run_dirs[0])[0])[0] + '/*.pdb')[0]
            print 'reference pdb: %s' % random_pdb
            mol = sasmol.SasMol(0)
            mol.read_pdb(random_pdb)

            if do_align:
                align_inputs = align.inputs()
                ref = sasmol.SasMol(0)
                ref.read_pdb(random_pdb)
                align_inputs.aa_goal = ref
                align_basis = ('((name[i] == "CA") and (segname[i] == "3H2A") '
                               'and (resid[i] > 5) and (resid[i] < 115))')
                align_inputs.goal_basis = align_basis
                align_inputs.move_basis = align_basis

            if cutoff:
                # keep structures below the cutoff value
                dcd_name = 'best_lt%s%s.dcd' % (str(cutoff), label)
            else:
                # keep the first N structures
                dcd_name = 'best_%d%s.dcd' % (N, label)
                result_df.sort(key, inplace=True)
                cutoff = result_df.iloc[N-1:N+1][key].mean()
            best_df = result_df.loc[result_df[key] < cutoff]
            iq_dir = op.join(data_file, 'foxs/')
            dcd_dir = op.join(data_file, 'monte_carlo/')
            mkdir_p(iq_dir)
            mkdir_p(dcd_dir)
            dcd_out = mol.open_dcd_write(op.join(dcd_dir, dcd_name))

            for (i, run_dir) in enumerate(run_dirs):
                # create read dcd pointer
                dcd_filename = glob.glob(run_dir + 'monte_carlo/*.dcd')
                assert len(dcd_filename) == 1, ('ERROR: %d dcd files in "%s"' %
                                                (len(dcd_filename), run_dir))
                dcd_in = mol.open_dcd_read(dcd_filename[0])

                for j in xrange(len(df_list[i])):
                    this_rg = df_list[i]['Rg'].iloc[j]
                    accept = df_list[i][key].iloc[j] < cutoff
                    index += 1
                    txt_file.write('%d\t%0.6f\t%0.6f\n' %
                                   (index, this_rg, accept))
                    # read specified dcd frame then save to dcd output
                    mol.read_dcd_step(dcd_in, df_list[i]['id'].iloc[j])
                    if accept:
                        n_accept += 1
                        if do_align:
                            # aligning here is much slower than
                            # aligning all frames at once
                            align_inputs.aa_move = mol
                            align.align_mol(align_inputs)

                        mol.write_dcd_step(dcd_out, 0, n_accept)
                        # write the iq data to a folder
                        best_iq = data_iq_list[i][:, [0, j+1]]
                        fname = 'foxs_%05d.dat' % n_accept
                        np.savetxt(op.join(iq_dir, fname), best_iq)

                mol.close_dcd_read(dcd_in[0])
            mol.close_dcd_write(dcd_out)

        else:
            for (i, run_dir) in enumerate(run_dirs):
                for j in xrange(len(df_list[i])):
                    this_rg = df_list[i]['Rg'].iloc[j]
                    accept = df_list[i][key].iloc[j] < cutoff
                    index += 1
                    txt_file.write('%d\t%0.6f\t%0.6f\n' %
                                   (index, this_rg, accept))


def plot_discrepancy(x2rg_df, all_data_iq, goal_iq, data_file, prefix='',
                     residual = r'$\chi^2$', key='x2', sub_label=True, N=None):

    n_total = len(x2rg_df)
    n_best = max(int(n_total * 0.1), 3)
    x2rg_best = x2rg_df.sort(key)[:n_best]

    plt.figure(figsize=(9, 4.5))
    gs = GridSpec(8, 2, hspace=0)

    ax0 = plt.subplot(gs[:, 0])
    ax0.text(0.01, 0.01, '%d structures' % n_total, verticalalignment='bottom',
             horizontalalignment='left', transform=ax0.transAxes)
    if sub_label:
        ax0.text(-0.03, -0.15, '(a)', verticalalignment='bottom',
                 horizontalalignment='left', transform=ax0.transAxes)
        ax0.text(1.13, -0.15, '(b)', verticalalignment='bottom',
                 horizontalalignment='left', transform=ax0.transAxes)

    ax0.plot(x2rg_df['Rg'], x2rg_df[key], 'o', mec=gp.qual_color(0), mfc='none')
    # plt.xlim(rg_range)
    plt.ylabel(residual)
    plt.xlabel(r'$R_g\,(\AA)$')

    # get the best, worst and average I(Q)
    best_x2 = x2rg_df[key].min()
    best_series = x2rg_df[x2rg_df[key] == best_x2]
    i_best = best_series.index[0] + 1  # first column is the Q values
    if goal_iq[0, 0] < 0.00001:
        xlim_min = goal_iq[1, 0] * .9 # set reasonable xlim
    else:
        xlim_min = goal_iq[0, 0] * .9

    best = all_data_iq[:, i_best]
    worst_x2 = x2rg_df[key].max()
    worst_series = x2rg_df[x2rg_df[key] == worst_x2]
    i_worst = worst_series.index[0] + 1  # first column is the Q values
    worst = all_data_iq[:, i_worst]
    if not N:
        average = all_data_iq[:, 1:].mean(axis=1)
        average_label = r'Average of All Structures'
    else:
        # get the index for the best N structures
        tmp_df = x2rg_df.sort(key)
        cutoff = tmp_df.iloc[N-1:N+1][key].mean()
        best_N_df = tmp_df.loc[tmp_df[key] < cutoff]

        # average the best N structures
        average = all_data_iq[:,best_N_df.index+1].mean(axis=1) # +1 b/c 0 is Q

        average_label = r'Average of Best %d Structures' % N

    ax0.set_yscale('log')
    plt.axis('tight')


    ax_dummy = plt.subplot(gs[:, 1])


    ax1 = plt.subplot(gs[:-1, 1])
    ax2 = plt.subplot(gs[-1, 1])
    # plot errorbar in two parts to get label order correct
    ax1.plot(goal_iq[:, 0], goal_iq[:, 1], 'o', ms=8, mfc='none',
             mec=gp.qual_color(0), label='Experimental')
    ax1.errorbar(goal_iq[:, 0], goal_iq[:, 1], goal_iq[:, 2], fmt=None,
                 ecolor=gp.qual_color(0))
    ax2.plot(goal_iq[:, 0], goal_iq[:, 1]*0, '--', c=gp.qual_color(0))

    ax1.plot(all_data_iq[:, 0], worst[:], c=gp.qual_color(3), linewidth=2,
             label=(r'Worst Structure'))
    ax2.plot(all_data_iq[:, 0], goal_iq[:, 1]-worst[:],
             c=gp.qual_color(3), linewidth=2,
             label=(r'Worst Structure'))

    ax1.plot(all_data_iq[:, 0], average[:], c=gp.qual_color(2), linewidth=2,
             label=average_label)
    ax2.plot(all_data_iq[:, 0], goal_iq[:, 1]-average[:],
             c=gp.qual_color(2), linewidth=2,
             label=average_label)

    ax1.plot(all_data_iq[:, 0], best[:], c=gp.qual_color(1), linewidth=2,
             label='Best Structure')
    ax2.plot(all_data_iq[:, 0], goal_iq[:, 1]-best[:],
             c=gp.qual_color(1), linewidth=2,
             label='Best Structure')


    ax1.set_ylabel(r'$I(Q)$')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    plt.axis('tight')
    ax1.set_xlim([xlim_min, 0.21])
    xlim = ax1.get_xlim()
    lg = ax1.legend(loc=3, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])

    ax2.set_xlabel(r'$Q\,(\AA^{-1})$')
    ax2.xaxis.labelpad = -1.3
    ax2.get_yaxis().set_ticks([])
    ax2.set_xscale('log')
    gp.zoomout(ax0, 0.2)
    ax2.set_xlim(xlim)
    # lg = ax2.legend(scatterpoints=1, numpoints=1, bbox_to_anchor=(1,1))
    # ax2.text(-0.03, -1.5, '(b)', verticalalignment='bottom',
             # horizontalalignment='left', transform=ax2.transAxes)

    plt.tight_layout()
    show = False
    if show:
        plt.show()
    else:
        if not sub_label:
            data_file += '_noLabels'
        if prefix:
            data_file = '%s_%s' % (prefix, data_file)
        save_name = '%s_fit' % (data_file)
        fig_file_name = op.join(os.getcwd(), save_name)
        # plt.savefig(fig_file_name[:-3] + 'png')
        plt.savefig(fig_file_name + '.png', dpi=400, bbox_inches='tight')
        if n_total > 20000:
            print 'View fit plot: \neog %s.png &' % fig_file_name
        else:
            print 'View fit plot: \nevince %s.eps &' % fig_file_name
            plt.savefig(fig_file_name + '.eps', dpi=400, bbox_inches='tight')
        # plt.show()
        plt.close()
    # plot_x2_components(goal_iq, all_data_iq[:, [0, i_best]], show=show,
                        # prefix=(prefix + '_' + data_file), s=s, o=o)



def plot_run_best(x2rg_df, all_data_iq, goal_iq, data_file, prefix='',
                  residual = r'$\chi^2$', key='x2'):

    n_total = len(x2rg_df)
    n_best = max(int(n_total * 0.1), 3)
    x2rg_best = x2rg_df.sort(key)[:n_best]

    plt.figure(figsize=(9, 7))
    plt.suptitle(data_file, fontsize=14)

    ax = plt.subplot(221)
    plt.title('all %d structures' % n_total)
    ax.plot(x2rg_df['Rg'], x2rg_df[key], 'o', mec=gp.qual_color(0), mfc='none')
    # plt.xlim(rg_range)
    plt.ylabel(residual)
    plt.xlabel(r'$R_g$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())

    # get the best, worst and average I(Q)
    best_x2 = x2rg_df[key].min()
    best_series = x2rg_df[x2rg_df[key] == best_x2]
    i_best = best_series.index[0] + 1  # first column is the Q values
    if goal_iq[0, 0] < 0.00001:
        xlim_min = goal_iq[1, 0] * .9
    else:
        xlim_min = goal_iq[0, 0] * .9

    best = all_data_iq[:, i_best]
    worst_x2 = x2rg_df[key].max()
    worst_series = x2rg_df[x2rg_df[key] == worst_x2]
    i_worst = worst_series.index[0] + 1  # first column is the Q values
    worst = all_data_iq[:, i_worst]
    average = all_data_iq[:, 1:].mean(axis=1)
    ax.set_yscale('log')

    ax = plt.subplot(222)
    plt.title(r'best %s=%0.1f, worst %s=%0.1f' % (residual, best_x2,
                                                  residual, worst_x2))
    ax.errorbar(goal_iq[:, 0], goal_iq[:, 1], goal_iq[:, 2], fmt='o',
                label='exp', ms=8, mfc='none', c=gp.qual_color(0),
                mec=gp.qual_color(0))
    # ax.plot(all_data_iq[:,0], best[:], '-->', mfc='none', ms=8,
    ax.plot(all_data_iq[:, 0], best[:], '-', mfc='none', ms=8,
            c=gp.qual_color(1), mec=gp.qual_color(1), linewidth=2,
            label='best (%d)' % i_best)
    # ax.plot(all_data_iq[:,0], average[:], '-.s', mfc='none', ms=8,
    ax.plot(all_data_iq[:, 0], average[:], '-', mfc='none', ms=8,
            c=gp.qual_color(2), mec=gp.qual_color(2), linewidth=2,
            label='average')
    # ax.plot(all_data_iq[:,0], worst[:], '-^', mfc='none', ms=8,
    ax.plot(all_data_iq[:, 0], worst[:], '-', mfc='none', ms=8,
            c=gp.qual_color(3), mec=gp.qual_color(3), linewidth=2,
            label='worst (%d)' % i_worst)
    plt.xlabel(r'$Q (\AA^{-1})$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    plt.ylabel(r'$I(Q)$')
    plt.yscale('log')
    plt.xscale('log')
    plt.axis('tight')
    gp.zoomout(ax, 0.2)
    plt.xlim([xlim_min, 0.21])
    lg = plt.legend(loc=3, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)

    best_colors = [gp.qual_color(1), gp.qual_color(8), gp.qual_color(9)]

    ax = plt.subplot(223)
    plt.title('best %d structures' % (n_best))
    plt.plot(x2rg_best['Rg'], x2rg_best[key], 'o', mec=gp.qual_color(0)
             , mfc='none')
    # plt.plot(x2rg_best['Rg'], x2rg_best[key], '.')
    plt.plot(x2rg_best.iloc[0]['Rg'], x2rg_best.iloc[0][key], '>',
             mec=gp.qual_color(1), mfc=best_colors[0],
             markersize=8)
    plt.plot(x2rg_best.iloc[1]['Rg'], x2rg_best.iloc[1][key], 's',
             mec=gp.qual_color(2), mfc=best_colors[1],
             markersize=8)
    plt.plot(x2rg_best.iloc[2]['Rg'], x2rg_best.iloc[2][key], '^',
             mec=gp.qual_color(3), mfc=best_colors[2],
             markersize=8)
    plt.ylabel(residual)
    plt.xlabel(r'$R_g$')
    # ax.set_yscale('log')
    plt.axis('tight')
    gp.zoomout(ax, 0.1)

    # update the worst and average I(Q)
    # worst_x2 = x2rg_best.x2.max()
    # worst_series = x2rg_best[x2rg_best.x2 == worst_x2]
    # i_worst = worst_series.index[0] + 1 # first column is the Q values
    # worst = all_data_iq[:,i_worst]
    average = all_data_iq[:, 1:].mean(axis=1)
    i_1st = x2rg_best.index[0] + 1  # first column is the Q values
    i_2nd = x2rg_best.index[1] + 1  # first column is the Q values
    i_3rd = x2rg_best.index[2] + 1  # first column is the Q values
    assert i_1st == i_best, 'incorrectly indexing'

    ax = plt.subplot(224)
    plt.title(r'best 3 %s values = %0.1f, %0.1f, %0.1f' % (
        residual, best_x2, x2rg_best[key].iloc[1], x2rg_best[key].iloc[2]))
    plt.errorbar(goal_iq[:, 0], goal_iq[:, 1], goal_iq[:, 2], fmt='o',
                 label='exp', ms=8, mfc='none', c=gp.qual_color(0),
                 mec=gp.qual_color(0))
    # plt.plot(all_data_iq[:,0], average[:], '-.s', label='average')
    plt.plot(all_data_iq[:, 0], best[:], '-', c=best_colors[0], linewidth=2,
             label=r'$1^{st}$ (%d)' % i_best)
    plt.plot(all_data_iq[:, 0], all_data_iq[:, i_2nd], '-', c=best_colors[1],
             linewidth=2, label=r'$2^{nd}$ (%d)' % i_2nd)
    plt.plot(all_data_iq[:, 0], all_data_iq[:, i_3rd], '-', c=best_colors[2],
             linewidth=2, label=r'$3^{rd}$ (%d)' % i_3rd)
    # plt.plot(all_data_iq[:,0], best[:], '-->',
    # label=r'$1^{st}$ (%d)' % i_best)
    # plt.plot(all_data_iq[:,0], all_data_iq[:,i_2nd], '-s',
    # label=r'$2^{nd}$ (%d)' % i_2nd)
    # plt.plot(all_data_iq[:,0], all_data_iq[:,i_3rd], '-^',
    # label=r'$3^{rd}$ (%d)' % i_3rd)
    plt.xlabel(r'$Q (\AA^{-1})$')
    plt.ylabel(r'$I(Q)$')
    plt.yscale('log')
    plt.xscale('log')
    plt.axis('tight')
    gp.zoomout(ax, 0.1)
    # plt.xlim([-0.01, 0.21])
    lg = plt.legend(loc='upper right', scatterpoints=1, numpoints=1)
    lg.draw_frame(False)

    plt.tight_layout()
    show = False
    if show:
        plt.show()
    else:
        fig_file_name = op.join(os.getcwd(), '%s_%s_fit' % (prefix, data_file))
        # plt.savefig(fig_file_name[:-3] + 'png')
        plt.savefig(fig_file_name + '.png', dpi=400, bbox_inches='tight')
        if n_total > 20000:
            print 'View fit plot: \neog %s.png &' % fig_file_name
        else:
            print 'View fit plot: \nevince %s.eps &' % fig_file_name
            plt.savefig(fig_file_name + '.eps', dpi=400, bbox_inches='tight')
        # plt.show()
        plt.close()
    # plot_x2_components(goal_iq, all_data_iq[:, [0, i_best]], show=show,
                       # prefix=(prefix + '_' + data_file), s=s, o=o)


def plot_x2_components(rf_data, mt_data, prefix=None, show=False,
                       residual = r'$\chi^2$'):
    '''
    plot the point-by-point components of the X^2 summation
    '''
    if residual == r'$\chi^2$':
        x2, components = get_x2_components(rf_data, mt_data)
        components *= 100 / x2
    else:
        residual = r'$R$'
        x2, components = get_r_components(rf_data, mt_data)
        components *= 100 / x2

    colors = gp.color_order  # gp.qual_color

    plt.figure()
    if prefix:
        plt.suptitle(prefix)
    gs1 = GridSpec(4, 2, hspace=0, wspace=0)#, left=0.1, right=0.9, bottom=0.075, top=0.925,

    # # log-log scale # #
    ax1 = plt.subplot(gs1[:3, 0])
    ax1.plot(rf_data[:, 0], rf_data[:, 1], 'o', ms=8, mfc='none',
             mec=colors(0), label='exp')
    # label='exp', ms=8, mfc='none', color='b')
    ax1.errorbar(rf_data[:, 0], rf_data[:, 1], rf_data[:, 2], fmt="none",
                 ecolor=colors(0))
    ax1.plot(mt_data[:, 0], mt_data[:, 1], '-', mfc='none', ms=8,
             c=colors(1), linewidth=2, label=(r'best %s= %0.1f' % (residual,
                                                                   x2)))
    plt.ylabel(r'$I(Q)$')
    plt.yscale('log')
    plt.xscale('log')
    plt.axis('tight')
    gp.zoomout(ax1, 0.1)
    ax1.get_yaxis().set_ticks([])
    ax1.get_xaxis().set_ticks([])
    # ax1.yaxis.labelpad = -.5
    lg = plt.legend(loc=3, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)

    ax2 = plt.subplot(gs1[3, 0])
    width_vals = np.copy(mt_data[:, 0])
    width_vals[:-1] = mt_data[1:, 0] - mt_data[:-1, 0]
    width_vals[-1] = width_vals[-2]
    rects = ax2.bar(
        mt_data[:, 0] - width_vals / 2, components, width=width_vals, color=colors(0))
    # for rect in rects:
        # height = rect.get_height()
        # label = int(round(height))
        # if label > 1:
            # ax2.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                     # str(label), ha='center', va='bottom')
    ax2.set_xscale('log')
    gp.zoomout(ax2, 0.7)
    ylim = ax2.get_ylim()
    ax2.set_ylim([0, ylim[1]])
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xlabel(r'$Q (\AA^{-1})$')
    ax2.set_ylabel(r'%% of %s' % residual)
    ax2.get_yaxis().set_ticks([])
    # ax2.xaxis.labelpad = -10
    # ax2.yaxis.labelpad = -.5

    # # log-lin scale # #
    ax3 = plt.subplot(gs1[:3, 1])
    ax3.plot(rf_data[:, 0], rf_data[:, 1], 'o', ms=8, mfc='none',
             mec=colors(0), label='exp')
    # label='exp', ms=8, mfc='none', color='b')
    ax3.errorbar(rf_data[:, 0], rf_data[:, 1], rf_data[:, 2], fmt="none",
                 ecolor=colors(0))
    ax3.plot(mt_data[:, 0], mt_data[:, 1], '-', mfc='none', ms=8,
             c=colors(1), linewidth=2, label=(r'best %s= %0.1f' % (residual,
                                                                   x2)))
    ax3.set_yscale('log')
    # ax3.set_xlim([-0.01, 0.21])
    ax3.set_ylim(ax1.get_ylim())
    ax3.get_xaxis().set_ticks([])
    ax3.get_yaxis().set_ticks([])
    # ax3.yaxis.labelpad = -.5

    ax4 = plt.subplot(gs1[3, 1])
    width_vals = np.copy(mt_data[:, 0])
    width_vals[:-1] = mt_data[1:, 0] - mt_data[:-1, 0]
    width_vals[-1] = width_vals[-2]
    rects = ax4.bar(
        mt_data[:, 0] - width_vals / 2, components, width=width_vals, color=colors(0))
    for rect in rects:
        height = rect.get_height()
        label = int(round(height))
        # if label > 1:
            # ax4.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                     # str(label), ha='center', va='bottom')
    # ax4.set_xscale('log')
    ax4.set_ylim(ax2.get_ylim())
    ax4.set_xlim(ax3.get_xlim())
    ax4.set_xlabel(r'$Q (\AA^{-1})$')
    ax4.get_yaxis().set_ticks([])
    # ax4.xaxis.labelpad = -10
    # ax4.yaxis.labelpad = -.5

    plt.tight_layout()
    if show:
        plt.show()
    else:
        out_file = 'x2_comp.eps'
        if prefix:
            out_file = '%s_%s' % (prefix, out_file)
        fig_file_name = op.join(os.getcwd(), out_file)
        print 'View x2 component plot: \nevince %s &' % fig_file_name
        plt.savefig(fig_file_name[:-3] + 'png', dpi=400, bbox_inches='tight')
        plt.savefig(fig_file_name, dpi=400, bbox_inches='tight')
        plt.close()

def auto_crop(img):
    if img.max() == 255:
        white = 255 * 3
    elif img.max() == 1.0:
        white = 1.0 * 3
    else:
        print 'WARNING: no white pixels in image, error possible'
        return img

    r, c = img.shape[0:2]
    val_white_row = white * c
    val_white_col = white * r
    white_rows = [i for i in xrange(r) if img[i, :].sum() == val_white_row]
    white_cols = [j for j in xrange(c) if img[:, j].sum() == val_white_col]

    return np.delete(np.delete(img, white_cols, axis=1), white_rows, axis=0)


def auto_crop_group(images):
    n_images = len(images)
    white_rows = [None] * n_images
    white_cols = [None] * n_images
    r = [None] * n_images
    c = [None] * n_images
    size = [None] * n_images

    for (k, img) in enumerate(images):
        if img.max() == 255:
            white = 255 * 3
        elif img.max() == 1.0:
            white = 1.0 * 3
        else:
            print 'WARNING: no white pixels in image %d, error possible', k
        r[k], c[k] = img.shape[0:2]
        val_white_row = white * c[k]
        val_white_col = white * r[k]
        white_rows[k] = [i for i in xrange(r[k])
                         if img[i, :].sum() == val_white_row]
        white_cols[k] = [i for i in xrange(c[k])
                         if img[:, i].sum() == val_white_col]

    crop_rows = [i for i in white_rows[0] for wr in white_rows[1:] if i in wr]
    crop_cols = [i for i in white_cols[0] for wc in white_cols[1:] if i in wc]

    for (k, img) in enumerate(images):
        images[k] = np.delete(np.delete(img, crop_cols, axis=1),
                              crop_rows[k], axis=0)

    return images


def pub_plot(x2rg_df, all_data_iq, goal_iq, density_plots, inset_files=[],
             inset_loc=[], prefix='', i0=False, cutoff=None, show=False):
    n_total = len(x2rg_df)
    n_best = max(int(n_total * 0.1), 3)
    x2rg_best = x2rg_df.sort('x2')[:n_best]
    plt.close()
    colors = gp.color_order  # gp.qual_color

    inset_fontsize = 11
    default_fontsize = 13

    fig = plt.figure(figsize=(10, 3.5))
    # # ~~ code for adding and positioning a suptitle ~~
    # plt.suptitle(data_file, fontsize=14)
    # plt.subplots_adjust(left=0.125, right = 0.9, bottom = 0.1, top = 0.9,
    # wspace = 0.2, hspace = 0.2)

    gs1 = GridSpec(1, 2, left=0.075, right=0.75, wspace=0.1, hspace=0,
                   top=0.95)
    ax1 = plt.subplot(gs1[:, 0])
    best_wrst_titles = [r'best $\chi^2$ model', r'worst $\chi^2$ model']
    inset_images = auto_crop_group([plt.imread(inset_files[0]),
                                    plt.imread(inset_files[1])])
    for (i, img) in enumerate(inset_images):
        ax = plt_inset.add_inset(ax1, inset_loc[i], axisbg='None')
        mask = np.tile(np.atleast_3d(np.any(img != 255, axis=2)),
                       (1, 1, img.shape[2]))  # this should mask the white
        img = np.ma.masked_where(mask, img)
        ax.imshow(img, interpolation='none')
        ax.axis('off')
        ax.set_title(best_wrst_titles[i], fontsize=inset_fontsize, y=0.95)
        ax.patch.set_visible(False)  # hide the 'canvas'

    ax1.text(0.01, 0.01, '%d structures' % n_total, verticalalignment='bottom',
             horizontalalignment='left', transform=ax1.transAxes,
             fontsize=inset_fontsize)
    ax1.text(-0.03, -0.15, '(a)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax1.transAxes,
             fontsize=default_fontsize)
    ax1.plot(x2rg_df['Rg'], x2rg_df['x2'], 'o', mec=colors(0), mfc='none')

    # rg_range = [np.floor(x2rg_df['Rg'].min()), np.ceil(x2rg_df['Rg'].max())]
    # plt.xlim(rg_range)
    ax1.set_ylabel(r'$\chi^2$', fontsize=default_fontsize)
    ax1.set_xlabel(r'$R_g$', fontsize=default_fontsize)
    ax1.set_zorder(ax.get_zorder() + 1)  # put ax1 in front of ax
    ax1.patch.set_visible(False)  # hide the 'canvas'

    if i0:
        all_data_iq = all_data_iq[1:]
        goal_iq = goal_iq[1:]

    # get the best, worst and average I(Q)
    ax2 = plt.subplot(gs1[:, 1])
    ax2.text(-0.03, -0.15, '(b)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax2.transAxes,
             fontsize=default_fontsize)

    # plot errorbar in two parts to get label order correct
    ax2.plot(goal_iq[:, 0], goal_iq[:, 1], 'o', ms=8, mfc='none',
             mec=colors(0), label='exp')
    # label='exp', ms=8, mfc='none', color='b')
    ax2.errorbar(goal_iq[:, 0], goal_iq[:, 1], goal_iq[:, 2], fmt="none",
                 ecolor=colors(0))

    best_x2 = x2rg_df.x2.min()
    if 3 == all_data_iq.shape[1]:
        best = all_data_iq[:, 1]
    else:
        best_series = x2rg_df[x2rg_df.x2 == best_x2]
        i_best = best_series.index[0] + 1  # first column is the Q values
        best = all_data_iq[:, i_best]
    ax2.plot(all_data_iq[:, 0], best[:], '-', mfc='none', ms=8,
             c=colors(1), linewidth=2, label=(r'best $\chi^2$= %0.1f' % best_x2))

    # average = all_data_iq[:,1:].mean(axis=1)
    # ax2.plot(all_data_iq[:,0], average[:], '-', mfc='none', ms=8,
    # c=colors(3), linewidth=2,
    # label='average')

    worst_x2 = x2rg_df.x2.max()
    if 3 == all_data_iq.shape[1]:
        worst = all_data_iq[:, 2]
    else:
        worst_series = x2rg_df[x2rg_df.x2 == worst_x2]
        i_worst = worst_series.index[0] + 1  # first column is the Q values
        worst = all_data_iq[:, i_worst]
    ax2.plot(all_data_iq[:, 0], worst[:], '-', mfc='none', ms=8, c=colors(2),
             linewidth=2, label=(r'worst $\chi^2$= %0.1f' % worst_x2))

    plt.xlabel(r'$Q (\AA^{-1})$', fontsize=default_fontsize)
    plt.ylabel(r'$I(Q)$', fontsize=default_fontsize)
    plt.yscale('log')
    plt.xscale('log')
    plt.axis('tight')
    gp.zoomout(ax2, 0.1)
    ax2.get_yaxis().set_ticks([])
    ax2.xaxis.labelpad = -1
    ax2.yaxis.labelpad = -1
    lg = plt.legend(loc=3, scatterpoints=1, numpoints=1,
                    prop={'size': default_fontsize})
    lg.draw_frame(False)

    best_colors = [colors(1), colors(8), colors(9)]

    gs2 = GridSpec(2, 1, left=0.76, bottom=0.01, right=0.99, top=0.99,
                   wspace=0, hspace=0)
    ax3 = plt.subplot(gs2[0, 0])
    ax3.text(0.01, 0.05, '(c)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax3.transAxes,
             fontsize=default_fontsize)
    ax4 = plt.subplot(gs2[1, 0])
    ax4.text(0.01, 0.09, '(d)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax4.transAxes,
             fontsize=default_fontsize)

    img1, img2 = auto_crop_group([plt.imread(density_plots[0]),
                                  plt.imread(density_plots[1])])
    ax3.imshow(img1)
    ax4.imshow(img2)

    ax3.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax3.axis('off')
    ax4.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax4.axis('off')

    if show:
        plt.show()
    else:
        out_file = 'result.eps'
        if prefix:
            out_file = '%s_%s' % (prefix, out_file)
        fig_file_name = op.join(os.getcwd(), out_file)
        print 'View pub plot: \nevince %s &' % fig_file_name
        plt.savefig(fig_file_name[:-3] + 'png', dpi=400, bbox_inches='tight')
        plt.savefig(fig_file_name, dpi=400, bbox_inches='tight')
        plt.close()

def method_plot(x2rg_df, all_data_iq, goal_iq, density_plots,  example_plots,
                pdb_file_name, dcd_file_names, sas_folders, all_density_plot=[],
                inset_files=[], inset_loc=[], prefix='', i0=False, cutoff=None,
                show=False):
    import sassie.calculate.convergence_test as convergence_test

    n_total = len(x2rg_df)
    n_best = max(int(n_total * 0.1), 3)
    x2rg_best = x2rg_df.sort('x2')[:n_best]
    plt.close()
    colors = gp.qual_color

    # inset_fontsize = 11
    # default_fontsize = 14

    fig = plt.figure(figsize=(12, 10))
    # # ~~ code for adding and positioning a suptitle ~~
    # plt.suptitle(data_file, fontsize=14)
    # plt.subplots_adjust(left=0.125, right = 0.9, bottom = 0.1, top = 0.9,
    # wspace = 0.2, hspace = 0.2)

    gs1 = GridSpec(4, 8, left=0.01, right=0.99, bottom=0.01, top=0.99,
                   wspace=0.3, hspace=0.3)
    ax1 = plt.subplot(gs1[:2, :4])
    best_wrst_titles = [r'Best $\chi^2$ Model', r'Worst $\chi^2$ Model']
    inset_images = auto_crop_group([plt.imread(inset_files[0]),
                                    plt.imread(inset_files[1])])
    for (i, img) in enumerate(inset_images):
        ax_c = plt_inset.add_inset(ax1, inset_loc[i], axisbg='None')
        mask = np.tile(np.atleast_3d(np.any(img != 255, axis=2)),
                       (1, 1, img.shape[2]))  # this should mask the white
        img = np.ma.masked_where(mask, img)
        ax_c.imshow(img, interpolation='none')
        ax_c.axis('off')
        ax_c.set_title(best_wrst_titles[i], y=0.95)
        ax_c.patch.set_visible(False)  # hide the 'canvas'

    ax1.text(0.01, 0.015, '(a)   %d Structures' % n_total,
             verticalalignment='bottom', horizontalalignment='left',
             transform=ax1.transAxes)
    ax1.plot(x2rg_df['Rg'], x2rg_df['x2'], 'o', mec=colors(0), mfc='none')
    ax1.set_ylabel(r'$\chi^2$')
    ax1.set_xlabel(r'$R_g$')
    ax1.set_yscale('log')
    rg_range = [np.round(x2rg_df['Rg'].min()), np.round(x2rg_df['Rg'].max())]
    ax1.set_xlim(rg_range)
    # x2_range = [0.1, np.round(x2rg_df['x2'].max())]
    # ax1.set_ylim(x2_range)
    # ax1.xaxis.labelpad = -1
    ax1.set_zorder(ax_c.get_zorder() + 1)  # put ax1 in front of ax
    ax1.patch.set_visible(False)  # hide the 'canvas'

    if i0:
        all_data_iq = all_data_iq[1:]
        goal_iq = goal_iq[1:]

    # get the best, worst and average I(Q)
    ax2 = plt.subplot(gs1[:2, 4:])
    ax2.text(0.01, 0.015, '(b)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax2.transAxes)
             # fontsize=default_fontsize)

    # plot errorbar in two parts to get label order correct
    ax2.plot(goal_iq[:, 0], goal_iq[:, 1], 'o', ms=8, mfc='none',
             mec=colors(0), label='Experiment')
    # label='exp', ms=8, mfc='none', color='b')
    ax2.errorbar(goal_iq[:, 0], goal_iq[:, 1], goal_iq[:, 2], fmt="none",
                 ecolor=colors(0))

    best_x2 = x2rg_df.x2.min()
    if 3 == all_data_iq.shape[1]:
        best = all_data_iq[:, 1]
    else:
        best_series = x2rg_df[x2rg_df.x2 == best_x2]
        i_best = best_series.index[0] + 1  # first column is the Q values
        best = all_data_iq[:, i_best]
    ax2.plot(all_data_iq[:, 0], best[:], c=colors(1), linewidth=2,
             label=(r'Best $\chi^2$= %0.1f' % best_x2))

    worst_x2 = x2rg_df.x2.max()
    if 3 == all_data_iq.shape[1]:
        worst = all_data_iq[:, 2]
    else:
        worst_series = x2rg_df[x2rg_df.x2 == worst_x2]
        i_worst = worst_series.index[0] + 1  # first column is the Q values
        worst = all_data_iq[:, i_worst]
    ax2.plot(all_data_iq[:, 0], worst[:], c=colors(3), linewidth=2,
             label=(r'Worst $\chi^2$= %0.1f' % worst_x2))

    plt.xlabel(r'$Q (\AA^{-1})$')
    plt.ylabel(r'$I(Q)$')
    plt.yscale('log')
    plt.xscale('log')
    plt.axis('tight')
    gp.zoomout(ax2, 0.1)
    ax2.get_yaxis().set_ticks([])
    ax2.xaxis.labelpad = -5
    ax2.yaxis.labelpad = -1
    lg = plt.legend(loc='center left', scatterpoints=1, numpoints=1)
                    # prop={'size': default_fontsize})
    lg.draw_frame(False)

    # convergence plots
    ax_c = plt.subplot(gs1[2:, :4])
    ax_c.text(0.01, 0.015, '(c)', verticalalignment='bottom',
            horizontalalignment='left', transform=ax_c.transAxes)
    if all_density_plot:
        inset_image = auto_crop(plt.imread(all_density_plot))
        ax_i = plt_inset.add_inset(ax_c, [0.4, 0.05, 0.6, 0.6], axisbg='None')
        mask = np.tile(np.atleast_3d(np.any(inset_image != 255, axis=2)),
                       (1, 1, inset_image.shape[2]))  # this should mask the white
        inset_image = np.ma.masked_where(mask, inset_image)
        ax_i.imshow(inset_image, interpolation='none')
        ax_i.axis('off')
        ax_i.set_title('All versus Best', y=1.05)
        ax_i.patch.set_visible(False)  # hide the 'canvas'

    # real-space convergence
    list_new_voxels = []
    list_occupied_voxels = []
    convergence_test.count_spatial_voxels(pdb_file_name, dcd_file_names,
                                          list_new_voxels, list_occupied_voxels)
    n_structures = sum([len(new_voxels) for new_voxels in list_new_voxels])
    occupied_voxels = np.zeros((n_structures, 2))
    occupied_voxels[:, 0] = np.arange(n_structures)
    for i in xrange(len(dcd_file_names)):
        rows = list_occupied_voxels[i][:, 0]
        occupied_voxels[rows, 1] = list_occupied_voxels[i][:, 1]
    ax_c.plot(occupied_voxels[:, 0], occupied_voxels[:, 1],
            label='Real-Space Convergence')
    ax_c.set_ylim((0, occupied_voxels[-1,1]))
    ax_c.set_xlabel('Number of Structures')
    # Make the y-axis label and tick labels match the line color.
    ax_c.set_ylabel('Number of Occupied Spatial Voxels', color=colors(0))
    for tl in ax_c.get_yticklabels():
        tl.set_color(colors(0))

    # reciprocal-space convergence
    ax_c2 = ax_c.twinx()
    iq_low = []
    iq_high = []
    iq_all = []
    list_new_grids = []
    list_occupied_grids = []
    n_q, n_spec = convergence_test.load_iq(sas_folders, 'dat', iq_low, iq_high,
                                           iq_all)
    convergence_test.count_sas_grids(sas_folders, iq_low, iq_high, iq_all, n_q,
                                     n_spec, list_new_grids,
                                     list_occupied_grids, 100)
    total_spec = n_spec.sum()
    occupied_grids = np.zeros((total_spec, len(sas_folders)+1))
    occupied_grids[:, 0] = np.arange(total_spec)
    for i in xrange(len(sas_folders)):
        rows = list_occupied_grids[i][:, 0]
        occupied_grids[rows, 1] = list_occupied_grids[i][:, 1]
    ax_c2.plot(occupied_grids[1:, 0], occupied_grids[1:, 1], color=colors(1),
               label='SAXS Convergence')
    ax_c2.set_ylim((0, occupied_grids[-1, 1]))
    # Make the y-axis label and tick labels match the line color.
    ax_c2.set_ylabel('Number of Occupied SAXS Grids', color=colors(1))
    for tl in ax_c2.get_yticklabels():
        tl.set_color(colors(1))

    ax_c.set_xlim((0, occupied_voxels[-1, 0]))
    ax_c.set_zorder(ax_i.get_zorder() + 1)  # put ax in front of ax_i
    ax_c.patch.set_visible(False)  # hide the 'canvas'

    # best fit w/ density plot
    gs2 = GridSpec(2, 4, left=0.57, right=1.0, bottom=0.0, top=0.47,
                   wspace=0.0, hspace=0.0)
    ax3 = plt.subplot(gs2[0, 0])
    ax3.text(0.01, 0.05, '(d)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax3.transAxes)
             # fontsize=default_fontsize)
    ax4 = plt.subplot(gs2[1, 0])
    ax4.text(0.01, 0.09, '(e)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax4.transAxes)
             # fontsize=default_fontsize)

    img1, img2 = auto_crop_group([plt.imread(density_plots[0]),
                                  plt.imread(density_plots[1])])
    ax3.imshow(img1)
    ax4.imshow(img2)

    # example 1
    ax5 = plt.subplot(gs2[0, 1])
    ax5.text(0.01, 0.05, '(f)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax5.transAxes)
    ax6 = plt.subplot(gs2[1, 1])
    ax6.text(0.01, 0.09, '(g)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax6.transAxes)

    ex1a, ex1b = auto_crop_group([plt.imread(example_plots[0]),
                                  plt.imread(example_plots[1])])
    ax5.imshow(ex1a)
    ax6.imshow(ex1b)

    # example 2
    ax7 = plt.subplot(gs2[0, 2])
    ax7.text(0.01, 0.05, '(h)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax7.transAxes)
    ax8 = plt.subplot(gs2[1, 2])
    ax8.text(0.01, 0.09, '(i)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax8.transAxes)

    ex2a, ex2b = auto_crop_group([plt.imread(example_plots[2]),
                                  plt.imread(example_plots[3])])
    ax7.imshow(ex2a)
    ax8.imshow(ex2b)

    # example 3
    ax9 = plt.subplot(gs2[0, 3])
    ax9.text(0.01, 0.05, '(j)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax9.transAxes)
    ax10 = plt.subplot(gs2[1, 3])
    ax10.text(0.01, 0.09, '(k)', verticalalignment='bottom',
             horizontalalignment='left', transform=ax10.transAxes)

    ex3a, ex3b = auto_crop_group([plt.imread(example_plots[4]),
                                  plt.imread(example_plots[5])])
    ax9.imshow(ex3a)
    ax10.imshow(ex3b)

    ax3.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax3.axis('off')
    ax4.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax4.axis('off')
    ax5.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax5.axis('off')
    ax6.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax6.axis('off')
    ax7.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax7.axis('off')
    ax8.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax8.axis('off')
    ax9.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax9.axis('off')
    ax10.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    ax10.axis('off')

    if show:
        plt.show()
    else:
        out_file = 'result.eps'
        if prefix:
            out_file = '%s_%s' % (prefix, out_file)
        fig_file_name = op.join(os.getcwd(), out_file)
        print 'View pub plot: \nevince %s &' % fig_file_name
        plt.savefig(fig_file_name[:-3] + 'png', dpi=400, bbox_inches='tight')
        plt.savefig(fig_file_name, dpi=400, bbox_inches='tight')
        # plt.show()
        plt.close()
        print 'pause'

def save_output(all_data_files, all_result_dfs, all_data_iqs, all_goal_iqs,
                sassie_run_dir, metric='i1_r', prefix=''):
    for (i, data_file) in enumerate(all_data_files):

        best_residual = all_result_dfs[i][metric].min()
        best_series = all_result_dfs[i][all_result_dfs[i][metric] == best_residual]
        i_best = best_series.index[0] + 1  # first column is the Q values

        wrst_residual = all_result_dfs[i][metric].max()
        wrst_series = all_result_dfs[i][all_result_dfs[i][metric] == wrst_residual]
        i_worst = wrst_series.index[0] + 1  # first column is the Q values

        all_result_dfs[i].index.name = 'index'
        all_result_dfs[i].to_csv(prefix + data_file + '_residulalRg.csv', sep=',')
        np.savetxt(prefix + data_file + '_bst_wrst.iq',
                   all_data_iqs[i][:, [0, i_best, i_worst]],
                   header='Q, best I(Q), worst I(Q)')
        np.savetxt(prefix + data_file + '_exp.iq', all_goal_iqs[i],
                   header='Q, I(Q), Error I(Q)')

        # save the best and worst dcd frame as a pdb
        mol = sasmol.SasMol(0)
        pdb_search = op.join(sassie_run_dir, op.split(
            best_series.iloc[0]['run'])[0], '*.pdb')
        random_pdb = glob.glob(pdb_search)[0]
        mol.read_pdb(random_pdb)

        best_mc_dir = op.join(sassie_run_dir, best_series.iloc[0]['run'],
                              'monte_carlo')
        best_dcd = glob.glob(op.join(best_mc_dir, '*.dcd'))
        assert len(best_dcd) == 1, 'ERROR: multiple dcd in %s' % best_mc_dir
        best_dcd = best_dcd[0]
        mol.read_single_dcd_step(best_dcd, int(best_series.iloc[0]['id']))
        mol.write_pdb(prefix + data_file + '_best.pdb', 0, 'w')

        wrst_mc_dir = op.join(sassie_run_dir, wrst_series.iloc[0]['run'],
                              'monte_carlo')
        wrst_dcd = glob.glob(op.join(wrst_mc_dir, '*.dcd'))
        assert len(wrst_dcd) == 1, 'ERROR: multiple dcd in %s' % wrst_mc_dir
        wrst_dcd = wrst_dcd[0]
        mol.read_single_dcd_step(wrst_dcd, wrst_series.iloc[0]['id'])
        mol.write_pdb(prefix + data_file + '_wrst.pdb', 0, 'w')

        # cp the psf-file to local directory
        random_psf = glob.glob(op.join(sassie_run_dir, op.split(
            best_series.iloc[0]['run'])[0], '*.psf'))[0]
        shutil.copy(random_psf, './')

if __name__ == '__main__':
    test = False
    if test:
        # import doctest
        # doctest.testmod()

        small = 0.3
        dummy_scale = 10
        dummy_offset = 0.5
        N = 50
        x = np.linspace(0.01, np.pi * 4, N)

        er = np.linspace(0.1, 0.5, N)
        y_rf = dummy_scale * np.sin(x) / x + dummy_offset
        rf_data = np.vstack([x, y_rf, er]).T

        y_in = np.sin(x) / x + (np.random.random(len(x)) - 0.5) * small
        in_data = np.vstack([x, y_in]).T

        compare_match(rf_data, in_data, dummy_scale, dummy_offset)

        y_rf = dummy_scale * np.sin(x) / x
        rf_data = np.vstack([x, y_rf, er]).T

        compare_match(rf_data, in_data, dummy_scale, 0)

    else:
        # data_rg_i0_df, data_files = examine_rg_i0(True)

        if False:
            # fig_rg_v_conc()
            # fig_sub_rg_v_conc(show=True)
            fig_sub_i0_v_conc(show=True)
            print 'wait'
            # fig_rg_v_salt(show=True)

        # tri5 = ['c500_3x167_k010',
                # 'c500_3x167_k050',
                # 'c500_3x167_k100',
                # 'c500_3x167_k200']
        # tet5 = ['c500_4x167_k010',
                # 'c500_4x167_k050',
                # 'c500_4x167_k100']
        # full
        tri_0 = ['c000_3x167_k010',
                 'c000_3x167_k050',
                 'c000_3x167_k100',
                 'c000_3x167_k200']
        tet_0 = ['c000_4x167_k010',
                'c000_4x167_k050',
                'c000_4x167_k100',
                'c000_4x167_mg1']
        di_0 = ['c000_2x167_k010']
        h5_0 = ['c000_4x167_h5_k010',
                'c000_4x167_h5_mg1']

        tri_all = tri_0 + [
            'c068_3x167_k010',
            'c125_3x167_k010',
            'c125_3x167_k050',
            'c125_3x167_k100',
            'c125_3x167_k200',
            'c250_3x167_k010',
            'c250_3x167_k050',
            'c250_3x167_k100',
            'c250_3x167_k200',
            'c500_3x167_k010',
            'c500_3x167_k050',
            'c500_3x167_k100',
            'c500_3x167_k200']
        # tri_all = tri_all[-4:]
        tet_all = tet_0 + []
        di_all = di_0 + []
        h5_all = h5_0 + [
            'c200_4x167_h5_k010',
            'c200_4x167_h5_mg1',
            'c400_4x167_h5_k010',
            'c400_4x167_h5_mg1',
            'c800_4x167_h5_k010',
            'c800_4x167_h5_mg1']

        tri_k010 = ['c000_3x167_k010']
        tri_k050 = ['c000_3x167_k050']
        tri_k100 = ['c000_3x167_k100']
        tri_k200 = ['c000_3x167_k200']

        data_dir = ('/home/schowell/data/Dropbox/gw_phd/paper_tetranucleosome/'
                    '1406data/chess/iqdata/coarse_grid/')
        data_ext = '_29.wi0'

        data_files = {'di': di_0, 'tri': tri_0, 'tet': tet_0, 'h5': h5_0}
        # data_files = {'di': di_all, 'tri': tri_all, 'tet': tet_all,
                      # 'h5': h5_all}

        # data_files['tri'] = tri_k010
        # maxx2 =  12.4301 # 3x167_k010
        # data_files['tri'] = tri_k050
        # maxx2 =  52.3015 # 3x167_k050
        # data_files['tri'] = tri_k100
        # maxx2 = 807.1556 # 3x167_k100
        # data_files['tri'] = tri_k200
        maxx2 =  12.7109 # 3x167_k200


        sassie_run_dir = '/home/schowell/data/myData/sassieRuns/'
        dimer_runs = []
        trimer_runs = []
        tetramer_runs = []
        dimer_runs += glob.glob(sassie_run_dir + '2x167*/run*/foxs/')

        k010_3x167_folders = glob.glob(sassie_run_dir + '3x167_k010/*/foxs/')
        k050_3x167_folders = glob.glob(sassie_run_dir + '3x167_k050/*/foxs/')
        k100_3x167_folders = glob.glob(sassie_run_dir + '3x167_k100/*/foxs/')
        mg01_3x167_folders = glob.glob(sassie_run_dir + '3x167_k200/*/foxs/')
        more_3x167_folders = glob.glob(sassie_run_dir + '3x167/*/foxs/')
        k010_3x167_folders.sort() # key=os.path.getmtime)
        k050_3x167_folders.sort() # key=os.path.getmtime)
        k100_3x167_folders.sort() # key=os.path.getmtime)
        mg01_3x167_folders.sort() # key=os.path.getmtime)
        more_3x167_folders.sort() # key=os.path.getmtime)
        trimer_runs = (k010_3x167_folders +
                       k050_3x167_folders +
                       k100_3x167_folders +
                       mg01_3x167_folders +
                       more_3x167_folders)

        k010_4x167_folders = glob.glob(sassie_run_dir + '4x167_k010/*/foxs/')
        k050_4x167_folders = glob.glob(sassie_run_dir + '4x167_k050/*/foxs/')
        k100_4x167_folders = glob.glob(sassie_run_dir + '4x167_k100/*/foxs/')
        mg01_4x167_folders = glob.glob(sassie_run_dir + '4x167_mg01/*/foxs/')
        more_4x167_folders = glob.glob(sassie_run_dir + '4x167/*/foxs/')
        k010_4x167_folders.sort() # key=os.path.getmtime)
        k050_4x167_folders.sort() # key=os.path.getmtime)
        k100_4x167_folders.sort() # key=os.path.getmtime)
        mg01_4x167_folders.sort() # key=os.path.getmtime)
        more_4x167_folders.sort() # key=os.path.getmtime)
        tetramer_runs = (k010_4x167_folders +
                         k050_4x167_folders +
                         k100_4x167_folders +
                         mg01_4x167_folders +
                         more_4x167_folders)

        dimer_runs = [run.replace('foxs/', '') for run in dimer_runs]
        trimer_runs = [run.replace('foxs/', '') for run in trimer_runs]
        tetramer_runs = [run.replace('foxs/', '') for run in tetramer_runs]

        for run in trimer_runs:
            # print run
            sub_dirs = ['sub%02d' % i for i in xrange(1, 100)]
            for sub_dir in sub_dirs:
                if op.exists(op.join(run, 'foxs/', sub_dir)):
                    trimer_runs.remove(run)
                    break

        run_dirs = {'di': dimer_runs, 'tri': trimer_runs, 'tet': tetramer_runs,
                    'h5': tetramer_runs}

        # array_types = ['di', 'tri', 'tet', 'h5']
        # array_types = ['di']
        # array_types = ['tet']
        array_types = ['tri']
        # array_types = ['h5']
        # array_types = ['h5', 'tri']

        best_dcd = False
        all_x2rg_dfs, i1_data_iqs, min_data_iqs, all_goal_iqs, all_data_files = (
            compare_iq(array_types, data_files, data_dir, data_ext, run_dirs,
                        cutoff=maxx2, best_dcd=best_dcd,
                        fresh=False))

        # best_dcd = True
        all_x2rg_dfs, i1_data_iqs, min_data_iqs, all_goal_iqs, all_data_files = (
            compare_iq(array_types, data_files, data_dir, data_ext, run_dirs,
                        cutoff=maxx2, best_dcd=best_dcd,
                        fresh=False))

        save_output(all_data_files, all_x2rg_dfs, i1_data_iqs, all_goal_iqs,
                    sassie_run_dir)

        best_dcd = False
        all_x2rg_dfs, i1_data_iqs, min_data_iqs, all_goal_iqs, all_data_files = (
            compare_iq(array_types, data_files, data_dir, data_ext, run_dirs,
                        cutoff=maxx2, best_dcd=best_dcd,
                        fresh=False))

        # best_dcd = True
        all_x2rg_dfs, i1_data_iqs, min_data_iq, all_goal_iqs, all_data_files = (
            compare_iq(array_types, data_files, data_dir, data_ext, run_dirs,
                        cutoff=maxx2, best_dcd=best_dcd,
                        fresh=False))

        save_output(all_data_files, all_x2rg_dfs, i1_data_iqs, all_goal_iqs,
                    sassie_run_dir)

        # density_plots = [['/home/schowell/data/code/pylib/x_dna/util/all/woI0ns50_s/3x167face_x2_lt_4p5_6A_voxels.png',
                          # '/home/schowell/data/code/pylib/x_dna/util/all/woI0ns50_s/3x167side_x2_lt_4p5_6A_voxels.png']]
        # best_wrst = [['/home/schowell/data/code/pylib/x_dna/util/all/woI0ns50_s/3x167_k200_best_face.tga',
                      # '/home/schowell/data/code/pylib/x_dna/util/all/woI0ns50_s/3x167_k200_wrst_side.tga']]
        # inset_loc = [[[0.1, 0.2, 0.3, 0.3], [0.5, 0.6, 0.3, 0.3]]]
        # pub_plot(all_x2rg_dfs[0], all_data_iqs[0], all_goal_iqs[0],
                 # density_plots[0], prefix=all_data_files[0],
                 # inset_files=best_wrst[0], inset_loc=inset_loc[0])

    print '\m/ >.< \m/'
