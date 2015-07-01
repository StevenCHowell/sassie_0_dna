#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Compare experimental data to the sassie structures
# Created: 20 March 2015
#
# $Id$
#
#000000001111111111222222222233333333334444444444555555555566666666667777777777
#234567890123456789012345678901234567890123456789012345678901234567890123456789


import logging
LOGGER = logging.getLogger(__name__) #add module name manually

import sys, os, glob, locale, errno, pickle, re
import os.path as op
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from scipy import interpolate
import matplotlib.pyplot as plt
import x_dna.util.gw_plot as gp
import matplotlib.gridspec as gridspec
import sassie.sasmol.sasmol as sasmol
import subprocess
debug = True

class MainError(Exception):
    pass

def pr():
    NotImplemented

def load_crysol(saspath, i0):
    '''
    load in the crysol output (Rg and I(Q))
    taken from "sassie/analyze/best.py"
    '''

    sasfilestring=saspath.split('/')[-2] + '_'

    nfilestring='find '+saspath+' -name "*.log" | grep -c log'
    nfout=os.popen(nfilestring,'r').readlines()
    nf=locale.atoi(nfout[0])
    print "\n# READING Rg VALUES FROM LOG FILES"
    log=[]

    for i in xrange(nf):

        mst = str(i+1).zfill(5)  #99999 maximum number of frames
        log.append(saspath+'/'+sasfilestring+mst)

    Rg = crysol_rg(nf, log)
    iq_data = load_crysol_iq(nf, log, i0)

    return Rg, iq_data, log

def load_foxs(saspath, i0=None):
    '''
    load in the FoXS output (Rg and I(Q))
    '''
    result_file = op.join(saspath, 'rg.csv')

    syn_files = glob.glob(op.join(saspath, '*.dat'))
    syn_files.sort()

    rg = [None]*len(syn_files)
    if op.isfile(result_file):
        rg_df = pd.read_csv(result_file, sep='\t')
        rg_df.index = rg_df['labels']
        save_rg = False
        for (i, syn_file) in enumerate(syn_files):
            rg[i] = rg_df['rg'].loc[op.split(syn_file)[1].split('.')[0]]
    else:
        try:
            dcd_file = glob.glob(op.join(op.join(op.split(saspath)[0],
                                                 'monte_carlo'),'*.dcd'))
            assert len(dcd_file) == 1, 'ERROR: not clear which dcd file to use'
            dcd_file = dcd_file[0]
            pdb_file = glob.glob(op.join(op.split(op.split(saspath)[0])[0],
                                         op.split(dcd_file)[1][:-4] + '*.pdb'))
            assert len(pdb_file) == 1, 'ERROR: not clear which pdb file to use'
            pdb_file = pdb_file[0]
            rg = dcd_rg(pdb_file, dcd_file)
        except:
            print 'did not calculate Rg using dcd and pdb'
        save_rg = True

    all_data = []
    iq = []
    q = []
    label = []
    for (i, syn_file) in enumerate(syn_files):
        data = np.loadtxt(syn_file)
        if i0:
            data[:,1:] /= data[0,1] * i0
        all_data.append(data)
        iq.append(data[:,1])
        q.append(data[:,0])
        label.append(op.basename(syn_file).split('.')[0])
        if not rg[i]:
            pdb_file = syn_file.replace('foxs', 'pdb').replace('.dat','')
            rg[i] = pdb_rg(pdb_file)

    if save_rg:
        rg_dict = {'rg': rg}
        rg_df = pd.DataFrame(rg_dict, index=label)
        rg_df.index.name = 'labels'
        rg_df.to_csv(result_file, sep='\t')
    rg = np.array(rg)

    iq_data = np.concatenate((q[0][..., None], np.array(iq).transpose()),axis=1)

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
    norm=1.0
    for i in xrange(ns):
        sinf=open(log[i]+'.int','r').readlines()
        nh=1
        sln=len(sinf)
        jj=0
        if 0 == i:
            # on the first iteration, setup the output array
            allspec=np.zeros((ns,(sln-nh)),np.float)
            qvar=np.zeros((sln-nh),np.float)
            print '# READING '+str(ns)+' SAS INT FILES'
        for j in xrange(sln):
            if i0:
                slin=sinf[j].split()
                if(j>=nh):
                    qval=locale.atof(slin[0])
                    qvar[jj]=qval
                    val1=locale.atof(slin[1])
                    if(jj==0):
                        norm=val1
                    nval1=i0*val1/norm
                    allspec[i][jj]=nval1
                    jj=jj+1
            else:
                slin=sinf[j].split()
                if(j>=nh):
                    qvar[jj]=locale.atof(slin[0])
                    allspec[i][jj] = locale.atof(slin[1])
                    jj += 1
    # not sure pandas is the way to go
    # iq_data = DataFrame(allspec, columns=qvar, index=range(1,ns+1)).T
    # iq_data.columns.name = 'id' ; iq_data.index.name = 'Q'

    # numpy seem to be the best option the for storing the I(Q) data
    iq_data = np.zeros(((sln-nh),ns+1))
    iq_data[:,0] = qvar
    iq_data[:,1:] = allspec.T

    return iq_data

def gnom_rg(gnom_files):
    '''
    read the Rg in from the crysol output
    taken from "sassie/analyze/best.py"
    '''
    prefix = []
    rg_reci=[] ; rg_real=[] ; rger_real=[]
    i0_reci=[] ; i0_real=[] ; i0er_real=[]
    for gnom_out in gnom_files:
        inf=open(gnom_out, 'r').readlines()
        ln=len(inf)
        for k in xrange(ln):
            lin=inf[k].split()
            # print lin
            if(len(lin)>0):
                if(lin[0] == 'Reciprocal'):
                    rg_reci.append(locale.atof(lin[4], func=float))
                    i0_reci.append(locale.atof(lin[8], func=float))
                elif(lin[0] =='Real' and lin[1] == 'space:'):
                    rg_real.append(locale.atof(lin[4], func=float))
                    i0_real.append(locale.atof(lin[6], func=float))
                    rger_real.append(locale.atof(lin[9], func=float))
                    i0er_real.append(locale.atof(lin[11], func=float))
                    prefix.append(op.split(
                        op.splitext(gnom_out)[0])[-1])


    gnom_dict = {'Rg grp':rg_reci, 'Rg grl':rg_real, 'RgEr grl':rger_real,
                 'I0 grp':i0_reci, 'I0 grl':i0_real, 'I0Er grl':i0er_real}
    return DataFrame(gnom_dict, index=prefix)

def crysol_rg(ns, log):
    '''
    read the Rg in from the crysol output
    taken from "sassie/analyze/best.py"
    '''
    rgarray=[] ; keep=[] ; rgcl=[] ; rgch=[]
    for i in xrange(ns):
        inf=open(log[i]+'.log','r').readlines()
        ln=len(inf)
        for k in xrange(ln):
            lin=inf[k].split()
            if(len(lin)>0):
                if(lin[0]=='Electron'):
                    trg=lin[3]
                    ftrg=locale.atof(trg)
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
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and op.isdir(path):
            pass
        else: raise

def compare_run_to_iq(run_dir, goal, ns, filter_dir, out_file=None,
                      q_base=1):

    assert op.exists(run_dir), 'No such run directory: %s' % run_dir
    if op.exists(op.join(run_dir,'crysol')) and False:
        syn_data_dir = op.join(run_dir,'crysol')
        ext = '/*.int'
    elif op.exists(op.join(run_dir,'foxs')):
        syn_data_dir = op.join(run_dir,'foxs')
        ext = '/*.dat'
    else:
        assert False, 'failed to find the calculated scattering data'
    syn_files = glob.glob(syn_data_dir + ext)
    assert len(syn_files) > 0, 'no I(Q) calculations in %s' % syn_data_dir
    syn_files.sort()

    # get the Rg and I(Q) for each structure from the calculation rusults
    if ext[-3:] == 'int':
        Rg, data_iq, labels = load_crysol(syn_data_dir, goal[0,1])
    elif ext[-3:] == 'dat':
        Rg, data_iq, labels = load_foxs(syn_data_dir, goal[0,1])

    data_iq, goal_iq = new_q_grid(ns, data_iq, goal, q_base=q_base,
                                  q_max=0.2)

    # # get the Kratky from the I(Q)
    # data_kratky = kratky(data_iq)
    # goal_kratky = kratky(goal_iq)
    # print data_kratky

    # create the data frame

    # compare I(Q) and/or Kratky of the experimental vs synthetic structures
    nf = data_iq.shape[1]
    matc_iq = np.zeros(data_iq.shape)
    s = np.zeros(nf-1)
    o = np.zeros(nf-1)
    X2 = np.zeros(nf-1)

    for i in xrange(1, nf):
        j = i-1
        # compare_match(goal_iq, data_iq[:,[0,i]], log=True)
        # X2 = (diff2 / sigma2).sum() # non-reduced chi-square distribution
        # X2 = (diff2 / goal).sum() # Pearson's X2 test-statistic (crazy units)
        matc_iq[:,[0,i]], s[j], o[j], X2[j] = scale_offset(data_iq[:,[0,i]],
                                                           goal_iq)

    res_dict = {'Rg':Rg, 'X2':X2, 'scale':s, 'offset':o, 'labels':labels}
    result_df = DataFrame(res_dict, index=range(1,nf))
    result_df.index.name = 'id'
    iq_df = DataFrame(matc_iq, columns=['Q']+labels)
    # save output to filter directory
    mkdir_p(filter_dir)
    if not out_file:
        out_file = op.join(filter_dir, 'rg_x2.out')
    result_df.to_csv(out_file, float_format='%5.10f', sep='\t')
    np.savetxt(op.join(filter_dir, 'goal.iq'), goal_iq) # small file
    np.save(op.join(filter_dir, 'data_iq.npy'), matc_iq)
    # np.savetxt('data.iq', data_iq) # too big to be useful as text file

    return result_df, matc_iq, iq_df, goal_iq

def new_q_grid(ns, data_iq, goal, q_max=0.2, q_base=1):
    ns -= 1
    q_min = goal[1, 0]
    new_q = polyspace(q_min, q_max, q_base, ns)

    new_q = np.concatenate(([0], new_q))
    # interpolate calculated data to be on the intended grid
    if len(new_q) == len(data_iq[:,0]) and np.allclose(data_iq[:,0], new_q):
        new_q = data_iq[:,0]
    else:
        if new_q[0] < data_iq[0,0]:
            new_q[0] = data_iq[0,0]
            print 'changing q_min from %f to %f' % (new_q[0], data_iq[0, 0])
        if new_q[-1] > data_iq[-1,0]:
            print 'changing q_max from %f to %f' % (new_q[-1], data_iq[-1, 0])
            new_q[-1] = data_iq[-1,0]
        interp_c = interpolate.interp1d(data_iq[:,0], data_iq[:,1:], axis=0,
                                             kind='cubic')
        data_int = np.zeros((len(new_q), len(data_iq[0,:])))
        data_int[:,0] = new_q
        data_int[:,1:] = interp_c(new_q)
        data_iq = data_int

    # interpolate reference data to be on the same grid as the crysol I(Q)
    interp_iq = interpolate.interp1d(goal[:,0], goal[:,1:], kind='cubic',
                                     axis=0)
    goal_iq = np.zeros((len(new_q), 3))
    goal_iq[:,0] = new_q
    goal_iq[:,1:] = interp_iq(new_q)
    return data_iq, goal_iq

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
    grid = np.linspace(x1 ** (1.0/p), x2 ** (1.0/p), n) ** p
    return grid

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
    assert (in_data[:,0]-rf_data[:,0]).sum() == 0, ('mismatch between input and'
                                                    ' reference x-grid')
    try:
        weights = 1/rf_data[:,2]
    except:
        weights = None
    offset, scale = np.polynomial.polynomial.polyfit(
        in_data[:,1], rf_data[:,1], 1, w=weights)
    mt_data = np.vstack([in_data[:,0], scale * in_data[:,1] + offset]).T

    X2 = get_X2(rf_data, mt_data)

    return mt_data, scale, offset, X2

def scale(in_data, rf_data):
    """
    determine the scale to match the input data to the reference
    data by minimizing the X2 calculation

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
    X2:      X^2 comparison between the reference data and matched input data

    See also
    --------
    match_poly, match_lstsq, scale_offset

    """
    assert (in_data[:,0]-rf_data[:,0]).sum() == 0, ('mismatch between input and'
                                                    ' reference x-grid')

    sigma2 = rf_data[:,2] * rf_data[:,2]
    scale = ( (rf_data[:,1] * in_data[:,1] / sigma2).sum() /
              (in_data[:,1] * in_data[:,1] / sigma2).sum() )

    mt_data = np.vstack([in_data[:,0], scale * in_data[:,1]]).T

    X2 = get_X2(rf_data, mt_data)

    return mt_data, scale, X2

def offset(in_data, rf_data):
    """
    determine the offset to match the input data to the reference
    data by minimizing the X2 calculation

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
    X2:      X^2 comparison between the reference data and matched input data

    See also
    --------
    match_poly, match_lstsq, scale_offset, scale

    """
    assert (in_data[:,0]-rf_data[:,0]).sum() == 0, ('mismatch between input and'
                                                    ' reference x-grid')

    sigma2 = rf_data[:,2] * rf_data[:,2]
    a = ( rf_data[:,1] / sigma2 ).sum()
    b = ( in_data[:,1] / sigma2 ).sum()
    c = ( 1 / sigma2 ).sum()
    offset = (a - b) / c

    mt_data = np.vstack([in_data[:,0], in_data[:,1] + offset]).T

    X2 = get_X2(rf_data, mt_data)

    return mt_data, offset, X2

def get_X2(rf_data, mt_data):
    diff = mt_data[:,1] - rf_data[:,1]
    diff2 = diff * diff
    er2 = rf_data[:,2] * rf_data[:,2]
    X2 = (diff2 / er2).sum()/len(rf_data)
    return X2

def scale_offset(in_data, rf_data):
    """
    determine the scale and offset to match the input data to the reference
    data by minimizing the X2 calculation
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
    X2:      X^2 comparison between the reference data and matched input data

    See also
    --------
    match_poly, match_lstsq, scale

    """
    e = 1E-4 # small parameter
    assert np.allclose((in_data[:,0]-rf_data[:,0]).sum(), 0, atol=e), (
        'mismatch between input and reference x-grid')

    sigma2 = rf_data[:,2] * rf_data[:,2]
    a = ( rf_data[:,1] / sigma2 ).sum()
    b = ( in_data[:,1] / sigma2 ).sum()
    c = ( 1 / sigma2 ).sum()
    d = ( rf_data[:,1] * in_data[:,1] / sigma2 ).sum()
    e = ( in_data[:,1] * in_data[:,1] / sigma2 ).sum()

    offset = (a*e - b*d)/(c*e - b*b)
    scale = (c*d - b*a)/(c*e - b*b)

    mt_data = np.vstack([in_data[:,0], scale * in_data[:,1] + offset]).T

    X2 = get_X2(rf_data, mt_data)

    return mt_data, scale, offset, X2

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
    assert (in_data[:,0]-rf_data[:,0]).sum() == 0, ('mismatch between input and'
                                                    ' reference x-grid')

    A = np.vstack([in_data[:,1], np.ones(len(in_data))]).T
    scale, offset = np.linalg.lstsq(A, rf_data[:,1])[0]
    mt_data = np.vstack([in_data[:,0], scale * in_data[:,1] + offset]).T

    X2 = get_X2(rf_data, mt_data)

    return mt_data, scale, offset, X2

def kratky(iq_data):
    '''
    iq_data should be a multi-column 2-D numpy array
    first column: Q
    returns a numpy array the same size where columns >= 2 are muliplied by q^2
    '''
    q = np.diag(iq_data[:,0])
    q2 = q.dot(q)
    kratky_data = np.zeros(iq_data.shape)
    kratky_data[:,0] = iq_data[:,0]
    kratky_data[:,1:] = q2.dot(iq_data[:,1:])

    return kratky_data

def compare_match(rf_data, in_data, dummy_scale=None, dummy_offset=None,
                  log=False):
    s = np.zeros(4)
    o = np.zeros(4)
    X2 = np.zeros(4)
    mt_data_polyf, s[0], o[0], X2[0] = match_poly(in_data, rf_data)
    mt_data_lstsq, s[1], o[1], X2[1] = match_lstsq(in_data, rf_data)
    mt_data_scale, s[2], X2[2]       = scale(in_data, rf_data)
    mt_data_sclof, s[3], o[3], X2[3] = scale_offset(in_data, rf_data)
    info = []
    for i in xrange(4):
        info.append('s=%0.1f, o=%0.3f, X2=%0.1f' % (s[i], o[i], X2[i]))
    if dummy_scale and dummy_offset:
        ref_str = ', s=%0.1f, o=%0.1f' % (dummy_scale, dummy_offset)
    else:
        ref_str = ''

    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.errorbar(rf_data[:,0], rf_data[:,1], yerr=rf_data[:,2], fmt='o',
                 label='reference' + ref_str)
    plt.plot(mt_data_polyf[:,0], mt_data_polyf[:,1], 's',
             label='poly, %s' %info[0])
    plt.plot(mt_data_sclof[:,0], mt_data_sclof[:,1], '^',
             label='sc-of, %s' %info[3])
    plt.plot(mt_data_lstsq[:,0], mt_data_lstsq[:,1], '<',
             label='lstsq, %s' %info[1])
    plt.plot(mt_data_scale[:,0], mt_data_scale[:,1], '>',
             label='scale, %s' %info[2])
    plt.plot(in_data[:,0], in_data[:,1], 'o', label='input')
    if log:
        plt.yscale('log')
        plt.xscale('log')
        plt.axis('tight')
        plt.legend(loc=3)
    else:
        plt.legend(loc=1)

    tl1 = plt.title("Comparison of different match methods")

    plt.subplot(2,1,2)
    ref = rf_data[:,1]
    plt.errorbar(rf_data[:,0], rf_data[:,1]-ref, yerr=rf_data[:,2], fmt='o',
                 label='reference')
    plt.plot(mt_data_polyf[:,0], mt_data_polyf[:,1]-ref, 's', label='poly')
    plt.plot(mt_data_sclof[:,0], mt_data_sclof[:,1]-ref, '^', label='sc-of')
    plt.plot(mt_data_lstsq[:,0], mt_data_lstsq[:,1]-ref, '<', label='lstsq')
    plt.plot(mt_data_scale[:,0], mt_data_scale[:,1]-ref, '>', label='scale')
    tl2 = plt.title('residuals')
    if log:
        plt.xscale('log')
        plt.axis('tight')
    plt.show()

    # i_sim = data_iq[1:,1]
    # i_exp = goal_iq[1:,1]
    # sigma = goal_iq[1:,2]
    # cs = np.linspace(0,20,201)
    # X2 = np.zeros(cs.shape)
    # dX2 = np.zeros(cs.shape)
    # for i, c in enumerate(cs):
        # X2[i] = ((c*i_sim-i_exp)**2/sigma**2).sum()/i_sim.shape[0]
        # dX2[i] = (2*I_sim*(c*i_sim-i_exp)/sigma**2).sum()/i_sim.shape[0]

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.hold()
    # plt.plot(cs, X2, label='X2')
    # plt.plot(cs, dX2, label='dX2/dc')
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
        plt.subplot(3,2,1)
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

        plt.subplot(3,2,3)
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

        plt.subplot(3,2,5)
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
        plt.subplot(3,2,2)
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

        plt.subplot(3,2,4)
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

        plt.subplot(3,2,6)
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

        plt.show()

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

    series.append(['N12merH5Mg1a_zeroCon', 'N12merH5Mg1ax4', 'N12merH5Mg1ax2', 'N12merH5Mg1ax1'])
    labels.append(r'12x167 H5 1mM $Mg^{2+}$')

    series.append(['N12merTE_zeroCon', 'N12merTEx4', 'N12merTEx2', 'N12merTEx1'])
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

    series.append(['N4merMg1_zeroCon', 'N4merMg1x4', 'N4merMg1x2', 'N4merMg1x1'])
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

    series.append(['tetraCtek010_zeroCon', 'tetraCtek010c012', 'tetraCtek010c050'])
    labels.append(r'4x167 10mM $K^+$ nc')

    series.append(['tetraCtek050_zeroCon', 'tetraCtek050c012', 'tetraCtek050c025',
                   'tetraCtek050c050'])
    labels.append(r'4x167 50mM $K^+$ nc')

    plt.figure()
    x_range = [-1, 1]

    for i in xrange(len(series)):
        plt.errorbar(df['conc'].loc[series[i]], df['Rg'].loc[series[i]],
                     df['RgEr'].loc[series[i]], label=labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i,'--'),
                     mec=gp.qual_color(i), mfc='none', ms=15)

    lg = plt.legend(loc=0, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'mg/mL')
    plt.title(r'$R_g$ vs mg/mL comparison')
    plt.xlim(x_range)
    plt.show()

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
    dod_labels.append(r'10mM $K^{+}$')

    dod.append(['c000_12x167_mg1',
                'c125_12x167_mg1',
                'c250_12x167_mg1',
                'c500_12x167_mg1'])
    dod_labels.append(r'1mM $Mg^{2+}$')

    dod.append(['c000_12x167_h5_k010',
                'c125_12x167_h5_k010',
                'c250_12x167_h5_k010',
                'c500_12x167_h5_k010'])
    dod_labels.append(r'H5 10mM $K^+$')

    dod.append(['c000_12x167_h5_mg1',
                'c125_12x167_h5_mg1',
                'c250_12x167_h5_mg1',
                'c500_12x167_h5_mg1'])
    dod_labels.append(r'H5 1mM $Mg^{2+}$')

    tri.append(['c000_3x167_k010',
                'c068_3x167_k010',
                'c125_3x167_k010',
                'c250_3x167_k010',
                'c500_3x167_k010'])
    tri_labels.append(r'10mM $K^+$')

    tri.append(['c000_3x167_k050',
                'c125_3x167_k050',
                'c250_3x167_k050',
                'c500_3x167_k050'])
    tri_labels.append(r'50mM $K^+$')
    inputs.run_name = 'run1'
    inputs.dcd_path = 'run1/dna_mc/'
    calc_foxs.main(inputs)
    print 'finished %s' % op.split(inputs.dcd_path[:-1])[0]

    inputs.run_name = 'run2'
    inputs.dcd_path = 'run2/dna_mc/'
    calc_foxs.main(inputs)
    print 'finished %s' % op.split(inputs.dcd_path[:-1])[0]

    tri.append(['c000_3x167_k100',
                'c125_3x167_k100',
                'c250_3x167_k100',
                'c500_3x167_k100'])
    tri_labels.append(r'100mM $K^+$')

    tri.append(['c000_3x167_k200',
                'c125_3x167_k200',
                'c250_3x167_k200',
                'c500_3x167_k200'])
    tri_labels.append(r'200mM $K^+$')

    di.append(['c000_2x167_k010',
               'c125_2x167_k010',
               'c250_2x167_k010',
               'c500_2x167_k010'])
    di_labels.append(r'10mM $K^+$')

    tet.append(['c000_4x167_k010',
                'c125_4x167_k010',
                'c250_4x167_k010',
                'c500_4x167_k010'])
    tet_labels.append(r'10mM $K^+$')

    tet.append(['c000_4x167_k050',
                'c125_4x167_k050',
                'c250_4x167_k050',
                'c500_4x167_k050'])
    tet_labels.append(r'50mM $K^+$')

    tet.append(['c000_4x167_k100',
                'c125_4x167_k100',
                'c250_4x167_k100',
                'c500_4x167_k100'])
    tet_labels.append(r'100mM $K^+$')

    tet.append(['c000_4x167_mg1',
                'c125_4x167_mg1',
                'c250_4x167_mg1',
                'c500_4x167_mg1'])
    tet_labels.append(r'1mM $Mg^{2+}$')

    tet.append(['c000_4x167_h5_k010',
                'c200_4x167_h5_k010',
                'c400_4x167_h5_k010',
                'c800_4x167_h5_k010'])
    tet_labels.append(r'H5 10mM $K^{+}$')

    tet.append(['c000_4x167_h5_mg1',
                'c200_4x167_h5_mg1',
                'c400_4x167_h5_mg1',
                'c800_4x167_h5_mg1'])
    tet_labels.append(r'H5 1mM $Mg^{2+}$')

    fig = plt.figure(figsize = (13, 10))
    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(hspace=0)
    x_range = [-0.05, 1.5]

    ## SUBPLOT(1,1) a)
    ax = plt.subplot(gs1[0])
    for i in xrange(len(tri)):
        plt.errorbar(df['conc'].loc[tri[i]], df['Rg'].loc[tri[i]],
                     df['RgEr'].loc[tri[i]], label=tri_labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i,'--'),
                     mec=gp.qual_color(i), mfc='none', ms=15, linewidth=2)

    lg = plt.legend(loc=4, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'mg/mL')
    # plt.title(r'3x167', x=0.3, y=0.92)
    ylim = np.array(plt.ylim())
    ylim[1] *= 1.08
    plt.ylim(ylim)
    tri_ylim = ylim # store for dimer ylim
    plt.xlim(x_range)
    ax.get_xaxis().set_ticks([])
    ax.text(0.03, 0.92, r'a) 3x167', verticalalignment='bottom', fontweight='bold',
        horizontalalignment='left', transform=ax.transAxes)

    ## SUBPLOT(2,1) b)
    ax = plt.subplot(gs1[2])
    for i in xrange(len(tet)):
        plt.errorbar(df['conc'].loc[tet[i]], df['Rg'].loc[tet[i]],
                     df['RgEr'].loc[tet[i]], label=tet_labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i,'--'),
                     mec=gp.qual_color(i), mfc='none', ms=15, linewidth=2)

    lg = plt.legend(loc=4, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'mg/mL')
    # plt.title(r'4x167', x=0.3, y=0.92)
    ylim = np.array(plt.ylim())
    ylim[1] *= 1.03
    plt.ylim(ylim)
    ax.set_yticks(ax.get_yticks()[:-2])
    plt.xlim(x_range)
    ax.text(0.03, 0.92, r'b) 4x167', verticalalignment='bottom', fontweight='bold',
        horizontalalignment='left', transform=ax.transAxes)


    x_range = [-0.05, 0.95]

    ## SUBPLOT(1,2) c)
    ax = plt.subplot(gs1[1])
    # plt.figure()
    for i in xrange(len(di)):
        this_color = gp.qual_color(i)
        plt.errorbar(df['conc'].loc[di[i]], df['Rg'].loc[di[i]],
                     df['RgEr'].loc[di[i]], label=di_labels[i],
                     c=this_color, fmt=gp.symbol_order(i,'--'),
                     mec=this_color, mfc='none', ms=15, linewidth=2)

    lg = plt.legend(loc=4, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    # plt.xlabel(r'mg/mL')
    # plt.title('2x167', x=0.3, y=0.92)
    # plt.text(.2, .9, r'2x167: $R_g$ vs mg/mL comparison',
             # horizontalalignment='center',nn
             # transform=plt.transAxes)
    plt.xlim(x_range)
    plt.ylim(tri_ylim)
    ax.get_xaxis().set_ticks([])
    ax.text(0.03, 0.92, r'c) 2x167', verticalalignment='bottom', fontweight='bold',
        horizontalalignment='left', transform=ax.transAxes)

    ## SUBPLOT(2,2) d)
    ax = plt.subplot(gs1[3])
    for i in xrange(len(dod)):
        plt.errorbar(df['conc'].loc[dod[i]], df['Rg'].loc[dod[i]],
                     df['RgEr'].loc[dod[i]], label=dod_labels[i],
                     c=gp.qual_color(i), fmt=gp.symbol_order(i,'--'),
                     mec=gp.qual_color(i), mfc='none', ms=15, linewidth=2)

    lg = plt.legend(loc=4, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'mg/mL')
    # plt.title(r'12x167', x=0.3, y=0.92)
    ylim = np.array(plt.ylim())
    ax.set_yticks(ax.get_yticks()[:-1])
    ylim[1] *= 1.03
    plt.ylim(ylim)
    plt.xlim(x_range)
    ax.text(0.03, 0.92, r'd) 12x167', verticalalignment='bottom', fontweight='bold',
            horizontalalignment='left', transform=ax.transAxes)

    fig.tight_layout()
    fig.savefig('Rg_v_mgmL.png')
    fig.savefig('Rg_v_mgmL.eps')
    if show: plt.show()

    return

def fig_rg_v_salt(show=False):
    df = load_rg_csv()
    tetA = ['tetraAtek010_zeroCon', 'tetraAtek050_zeroCon',
            'tetraAtek100_zeroCon']
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
    series = [di0, di5, tri0, tri5, tet0, tet5]
    labels = ['2x167 (0.0 mg/mL)',
              '2x167 (0.5 mg/mL)',
              '3x167 (0.0 mg/mL)',
              '3x167 (0.5 mg/mL)',
              '4x167 (0.0 mg/mL)',
              '4x167 (0.5 mg/mL)']
    fig = plt.figure()
    x_range = [-10, 220]
    for i in xrange(len(series)):
        plt.errorbar(df['KCl'].loc[series[i]], df['Rg'].loc[series[i]],
                     df['RgEr'].loc[series[i]], label=labels[i],
                     c = gp.qual_color(i), fmt=gp.symbol_order(i, '--'),
                     mec=gp.qual_color(i), mfc='none', ms=15)

    lg = plt.legend(loc=0, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)
    plt.ylabel(r'$R_g$')
    plt.xlabel(r'[$K^+$]')
    plt.title(r'$R_g$ vs [$K^+$] comparison')
    plt.xlim(x_range)
    if show: plt.show()
    fig.tight_layout()
    fig.savefig('Rg_v_salt.png')
    fig.savefig('Rg_v_salt.eps')
    return

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
    chess_matlab_file = 'rg_i0_chess_nor.csv'
    chess_matlab_df = pd.read_csv(op.join(chess_dir, chess_matlab_file))
    chess_matlab_df = chess_matlab_df.set_index(chess_matlab_df['prefix']).drop('prefix',1)

    chess_gnom_files = glob.glob(op.join(chess_dir, data_dir,'*.out'))
    chess_gnom_df = gnom_rg(chess_gnom_files)

    # nsls_dir = '/home/schowell/Dropbox/gw_phd/experiments/1406NSLS/'
    nsls_dir = '/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/1406data/nsls/'
    nsls_matlab_file = 'rg_i0_nsls_nor.csv'
    nsls_matlab_df = pd.read_csv(op.join(nsls_dir, nsls_matlab_file))
    nsls_matlab_df = nsls_matlab_df.set_index(nsls_matlab_df['prefix']).drop('prefix',1)

    nsls_gnom_files = glob.glob(op.join(nsls_dir, data_dir,'*.out'))
    nsls_gnom_df = gnom_rg(nsls_gnom_files)

    chess_df = pd.concat([chess_gnom_df, chess_matlab_df])
    nsls_df =  pd.concat([nsls_gnom_df, nsls_matlab_df])
    df = pd.concat([chess_df, nsls_df])
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

    return df

def evaluate_qqiq(array_types, data_files, data_dir, data_ext, run_dirs, ns):

    for array_type in array_types:
        for data_file in data_files[array_type]:
            full_file = op.join(data_dir, data_file + data_ext)
            assert op.exists(full_file), (
                    "No such file: '%s'" % full_file)
            data = np.loadtxt(full_file)

            df_list = [] ; data_qqiq_list = []
            for run_dir in run_dirs[array_type]:
                filter_dir = op.join(run_dir, data_file + '_filter')
                if op.exists(filter_dir + '/rg_x2.out'):
                    result_df = pd.read_csv(filter_dir + '/rg_x2.out', sep='\t')
                    data_qqiq = np.load(filter_dir + '/data_qqiq.npy')
                    goal_qqiq = np.loadtxt(filter_dir + '/goal.qqiq')
                else:
                    result_df, data_qqiq, goal_iq = compare_run_to_qqiq(
                        run_dir, data, ns[array_type], filter_dir)
                run_name = run_dir.split('/')[-3] + '/' + run_dir.split('/')[-2]
                result_df['run'] = run_name
                df_list.append(result_df)
                data_qqiq_list.append(data_qqiq)

            # combine result DataFrames and data_iq arrays
            x2rg_df = pd.concat(df_list)
            x2rg_df.index = range(len(x2rg_df))

            all_data_qqiq = data_qqiq_list[0][:,0:1]
            for data_qqiq in data_qqiq_list:
                all_data_qqiq = np.concatenate((all_data_qqiq, data_qqiq[:,1:]),
                                             axis=1)

            plot_run_best(x2rg_df, all_data_qqiq, goal_iq, data_file)

def evaluate_iq(array_types, data_files, data_dir, data_ext, run_dirs, ns,
                prefix='', do_plot='True', cutoff=None, join_dcd=False,
                q_base=1):
    all_x2rg_dfs = []
    all_iqs_dfs = []
    all_goal_iqs = []
    all_data_files = []

    if q_base:
        prefix += 'base%d_' % q_base
    else:
        prefix += 'lin_'

    for array_type in array_types:
        for data_file in data_files[array_type]:
            full_file = op.join(data_dir, data_file + data_ext)
            assert op.exists(full_file), (
                    "No such file: '%s'" % full_file)

            df_list = [] ; data_iq_list = []
            for run_dir in run_dirs[array_type]:
                filter_dir = op.join(run_dir, data_file + '_filter')
                if q_base:
                    out_file = op.join(filter_dir, 'rg_x2_base%d.out' %
                                       q_base)
                else:
                    out_file = op.join(filter_dir, 'rg_x2_lin.out')
                if op.exists(out_file) and False:
                    print 'loading rg and X^2 values from %s' % out_file
                    result_df = pd.read_csv(out_file, sep='\t')
                    data_iq = np.load(filter_dir + '/data_iq.npy')
                    goal_iq = np.loadtxt(filter_dir + '/goal.iq')
                else:
                    print 'loading rg then calculating X^2 values'
                    data = np.loadtxt(full_file)
                    result_df, data_iq, iq_df, goal_iq = compare_run_to_iq(
                        run_dir, data, ns[array_type], filter_dir, out_file,
                        q_base)
                run_name = run_dir.split('/')[-3] + '/' + run_dir.split('/')[-2]
                result_df['run'] = run_name
                df_list.append(result_df)
                data_iq_list.append(data_iq)
            if cutoff:
                write_filter_output(run_dirs[array_type], df_list, cutoff,
                                    join_dcd)

            # combine result DataFrames and data_iq arrays
            x2rg_df = pd.concat(df_list)
            x2rg_df.index = range(len(x2rg_df))
            x2rg_df.sort('X2', inplace=True)

            best5 = x2rg_df.iloc[:5]
            best5.to_csv(data_file + '_best.csv', float_format='%5.10f',
                         sep='\t')

            q = data_iq_list[0][:,:1]
            for (i, data_iq) in enumerate(data_iq_list):
                if i == 0:
                    all_data_iq = np.concatenate((q, data_iq[:,1:]), axis=1)
                else:
                    all_data_iq = np.concatenate((all_data_iq, data_iq[:,1:]),
                                                 axis=1)
            if do_plot:
                plot_run_best(x2rg_df, all_data_iq, goal_iq, data_file, prefix)

            all_x2rg_dfs.append(x2rg_df)
            all_goal_iqs.append(goal_iq)
            all_data_files.append(data_file)
            try:
                all_iqs_dfs.append(iq_df)
            except:
                print 'all_iqs_dfs contains no output'

    return  all_x2rg_dfs, all_iqs_dfs, all_goal_iqs, all_data_files

def write_filter_output(run_dirs, df_list, cutoff, join_dcd=False,
                        catdcd_exe='/home/myPrograms/bin/catdcd'):
    filter_dir = op.join(op.split(op.split(op.split(run_dirs[0])[0])[0])[0],
                         'filter')
    txt_name =  'rglowweights.txt'
    print 'saving output dcd and "%s" to %s' % (txt_name, filter_dir)
    mkdir_p(filter_dir)
    full_dcd_out = op.join(filter_dir,
                           'collected_%s.dcd' % filter_dir.split('/')[-2])
    tmp_dcd_out = op.join(filter_dir, 'tmp.dcd')
    index = 0
    with open(op.join(filter_dir, txt_name), 'w') as txt_file:
        txt_file.write('# file generated on FILL THIS IN\n')
        txt_file.write('# structure, Rg, weight\n')
        for (i, run_dir) in enumerate(run_dirs):
            for j in xrange(len(df_list[i])):
                this_rg = df_list[i]['Rg'].iloc[j]
                accept = df_list[i]['X2'].iloc[j] < cutoff
                index += 1
                txt_file.write('%d\t%0.6f\t%0.6f\n' % (index, this_rg, accept))
            if join_dcd:
                dcd_name = op.join(op.join(run_dir, 'foxs'), 'foxs_filtered.dcd')
                assert op.exists(dcd_name), ('ERROR!!!, could not find dcd file: %s'
                                             % dcd_file)
                if i == 0:
                    bash_cmd = 'cp %s %s' % (dcd_name, full_dcd_out)
                    subprocess.call(bash_cmd.split())
                else:
                    bash_cmd1 = ('%s -o %s %s %s' % (catdcd_exe, tmp_dcd_out,
                                                     full_dcd_out, dcd_name))
                    subprocess.call(bash_cmd1.split())
                    bash_cmd2 = 'mv %s %s' % (tmp_dcd_out, full_dcd_out)
                    subprocess.call(bash_cmd2.split())


def plot_run_best(x2rg_df, all_data_iq, goal_iq, data_file, prefix=''):

    n_total = len(x2rg_df)
    n_best = max(int(n_total * 0.1), 3)
    x2rg_best = x2rg_df.sort('X2')[:n_best]

    plt.figure(figsize=(9, 7))
    plt.suptitle(data_file, fontsize=14)

    # plt.subplots_adjust(left=0.125, right = 0.9, bottom = 0.1, top = 0.9,
                        # wspace = 0.2, hspace = 0.2)
    # plt.subplots_adjust(left=0.0, hspace=0.001)
    # rg_range = [np.floor(x2rg_df['Rg'].min()), np.ceil(x2rg_df['Rg'].max())]

    ax = plt.subplot(221)
    plt.title('all %d structures' % n_total)
    ax.plot(x2rg_df['Rg'], x2rg_df['X2'], 'o', mec='b', mfc='none')
    # plt.xlim(rg_range)
    plt.ylabel(r'$X^2$')
    plt.xlabel(r'$R_g$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())

    # get the best, worst and average I(Q)
    best_X2 = x2rg_df.X2.min()
    best_series = x2rg_df[x2rg_df.X2 == best_X2]
    i_best = best_series.index[0] + 1 # first column is the Q values
    best = all_data_iq[:,i_best]
    worst_X2 = x2rg_df.X2.max()
    worst_series = x2rg_df[x2rg_df.X2 == worst_X2]
    i_worst = worst_series.index[0] + 1 # first column is the Q values
    worst = all_data_iq[:,i_worst]
    average = all_data_iq[:,1:].mean(axis=1)

    ax = plt.subplot(222)
    plt.title(r'best $X^2$=%0.1f, worst $X^2$=%0.1f' % (best_X2, worst_X2))
    ax.errorbar(goal_iq[1:,0], goal_iq[1:,1], goal_iq[1:,2], fmt = '-o',
                label='exp', ms=8, mfc='none', c=gp.qual_color(0),
                mec=gp.qual_color(0))
    # ax.plot(all_data_iq[1:,0], best[1:], '-->', mfc='none', ms=8,
    ax.plot(all_data_iq[1:,0], best[1:], '-', mfc='none', ms=8,
            c=gp.qual_color(1), mec=gp.qual_color(1),
            label='best (%d)' % i_best)
    # ax.plot(all_data_iq[1:,0], average[1:], '-.s', mfc='none', ms=8,
    ax.plot(all_data_iq[1:,0], average[1:], '-', mfc='none', ms=8,
            c=gp.qual_color(2), mec=gp.qual_color(2),label='average')
    # ax.plot(all_data_iq[1:,0], worst[1:], '-^', mfc='none', ms=8,
    ax.plot(all_data_iq[1:,0], worst[1:], '-', mfc='none', ms=8,
            c=gp.qual_color(3), mec=gp.qual_color(3),
            label='worst (%d)' % i_worst)
    plt.xlabel(r'$Q (\AA^{-1})$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    plt.ylabel(r'$I(Q)$')
    plt.yscale('log')
    plt.xscale('log')
    plt.axis('tight')
    gp.zoomout(ax, 0.1)
    lg = plt.legend(loc=3, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)

    plt.subplot(223)
    plt.title('best %d structures' %(n_best))
    plt.plot(x2rg_best['Rg'], x2rg_best['X2'], 'o', mec='b', mfc='none')
    # plt.plot(x2rg_best['Rg'], x2rg_best['X2'], '.')
    plt.plot(x2rg_best.iloc[0]['Rg'], x2rg_best.iloc[0]['X2'], '>',
             mec=gp.qual_color(1), mfc=gp.qual_color(1),
             markersize=8)
    plt.plot(x2rg_best.iloc[1]['Rg'], x2rg_best.iloc[1]['X2'], 's',
             mec=gp.qual_color(2), mfc=gp.qual_color(2),
             markersize=8)
    plt.plot(x2rg_best.iloc[2]['Rg'], x2rg_best.iloc[2]['X2'], '^',
             mec=gp.qual_color(3), mfc=gp.qual_color(3),
             markersize=8)
    plt.ylabel(r'$X^2$')
    plt.xlabel(r'$R_g$')
    # plt.xlim(rg_range)

    # update the worst and average I(Q)
    # worst_X2 = x2rg_best.X2.max()
    # worst_series = x2rg_best[x2rg_best.X2 == worst_X2]
    # i_worst = worst_series.index[0] + 1 # first column is the Q values
    # worst = all_data_iq[:,i_worst]
    average = all_data_iq[:,1:].mean(axis=1)
    i_1st = x2rg_best.index[0] + 1 # first column is the Q values
    i_2nd = x2rg_best.index[1] + 1 # first column is the Q values
    i_3rd = x2rg_best.index[2] + 1 # first column is the Q values
    assert i_1st == i_best, 'incorrectly indexing'

    ax = plt.subplot(224)
    plt.title(r'best 3 $X^2$s = %0.1f, %0.1f, %0.1f' % (
        best_X2, x2rg_best.iloc[1].X2, x2rg_best.iloc[2].X2))
    plt.errorbar(goal_iq[1:,0], goal_iq[1:,1], goal_iq[1:,2], fmt = 'o',
                 label='exp', ms=8, mfc='none', c=gp.qual_color(0),
                mec=gp.qual_color(0))
    # plt.plot(all_data_iq[1:,0], average[1:], '-.s', label='average')
    plt.plot(all_data_iq[1:,0], best[1:], '-',
             label=r'$1^{st}$ (%d)' % i_best)
    plt.plot(all_data_iq[1:,0], all_data_iq[1:,i_2nd], '-',
             label=r'$2^{nd}$ (%d)' % i_2nd)
    plt.plot(all_data_iq[1:,0], all_data_iq[1:,i_3rd], '-',
             label=r'$3^{rd}$ (%d)' % i_3rd)
    # plt.plot(all_data_iq[1:,0], best[1:], '-->',
             # label=r'$1^{st}$ (%d)' % i_best)
    # plt.plot(all_data_iq[1:,0], all_data_iq[1:,i_2nd], '-s',
             # label=r'$2^{nd}$ (%d)' % i_2nd)
    # plt.plot(all_data_iq[1:,0], all_data_iq[1:,i_3rd], '-^',
             # label=r'$3^{rd}$ (%d)' % i_3rd)
    plt.xlabel(r'$Q (\AA^{-1})$')
    plt.ylabel(r'$I(Q)$')
    plt.yscale('log')
    plt.xscale('log')
    plt.axis('tight')
    gp.zoomout(ax, 0.1)
    lg = plt.legend(loc=3, scatterpoints=1, numpoints=1)
    lg.draw_frame(False)

    plt.tight_layout()
    # plt.show()

    fig_file_name = op.join(os.getcwd(), prefix + data_file + '_fit.eps')
    # plt.savefig(fig_file_name[:-3] + 'png')
    print 'storing fit plot as: %s' % fig_file_name
    plt.savefig(fig_file_name[:-3] + 'png', dpi=400, bbox_inches='tight')
    plt.savefig(fig_file_name, bbox_inches='tight')

    return


if __name__ == '__main__':
    test = False
    if test:
        # import doctest
        # doctest.testmod()

        small = 0.3; dummy_scale = 10; dummy_offset = 0.5; N = 50
        x = np.linspace(0.01, np.pi * 4, N)

        er = np.linspace(0.1, 0.5, N)
        y_rf = dummy_scale * np.sin(x)/x + dummy_offset
        rf_data = np.vstack([x, y_rf, er]).T

        y_in = np.sin(x)/x + (np.random.random(len(x)) - 0.5) * small
        in_data = np.vstack([x,y_in]).T

        compare_match(rf_data, in_data, dummy_scale, dummy_offset)

        y_rf = dummy_scale * np.sin(x)/x
        rf_data = np.vstack([x, y_rf, er]).T

        compare_match(rf_data, in_data, dummy_scale, 0)

    else:
        # fig_rg_v_conc()
        # data_rg_i0_df, data_files = examine_rg_i0(True)

        if False:
            fig_sub_rg_v_conc(show=True)
            fig_rg_v_salt(show=True)

        tri5 = ['c500_3x167_k010',
                'c500_3x167_k050',
                'c500_3x167_k100',
                'c500_3x167_k200']
        tet5 = ['c500_4x167_k010',
                'c500_4x167_k050',
                'c500_4x167_k100']
        # full
        tri0 = ['c000_3x167_k010',
                'c000_3x167_k050',
                'c000_3x167_k100',
                'c000_3x167_k200']
        tet0 = ['c000_4x167_k010',
                'c000_4x167_k050',
                'c000_4x167_k100',
                'c000_4x167_mg1']
        di0 = ['c000_2x167_k010']
        h5 = ['c000_4x167_h5_k010',
              'c000_4x167_h5_mg1']

        # quicker to run
        # tet0 = [tet0[-1]]
        # tet0 = ['c000_4x167_mg1']
        # tri0 = ['c000_3x167_k010']

        data_files = {'di': di0, 'tri': tri0, 'tet': tet0, 'h5': h5}

        sassie_run_dir = '/home/schowell/data/myData/sassieRuns/'
        dimer_runs = glob.glob(sassie_run_dir + 'dimer/flex*/run*/')
        trimer_runs = glob.glob(sassie_run_dir + 'trimer/flex*/run*/')
        tetramer_runs = glob.glob(sassie_run_dir + 'tetramer/flex*/run*/')
        # dimer_runs += glob.glob(sassie_run_dir + '2x167_*/run*/')
        # trimer_runs += glob.glob(sassie_run_dir + '3x167_*/run*/')
        tetramer_runs += glob.glob(sassie_run_dir + '4x167_*/run*/')
        run_dirs = {'di': dimer_runs, 'tri': trimer_runs, 'tet': tetramer_runs,
                    'h5': tetramer_runs}

        data_dir = ('/home/schowell/data/'
                    'Dropbox/gw_phd/paper_tetranucleosome/1406data/iqdata/')
        data_ext = '.i0q'

        array_types = ['di', 'tri', 'tet', 'h5']
        # array_types = ['tet']
        # array_types = ['tri']
        # array_types = ['h5']
        ns = {'di': 26, 'tri': 26,'tet': 26, 'h5': 26}  # N Q grid points
        evaluate_iq(array_types, data_files, data_dir, data_ext, run_dirs, ns,
                    cutoff=350, join_dcd=False, q_base=2)

    ############################ pseudo code ###########################
    # create a data frame containing the information for each structure
       # columns: Rg, X2 I(Q), X2 QI(Q), X2 P(r), Psi, Phi, h
       # rows: each structure

    # create a data frame of the crysol output
       # not sure if I should use columns or rows for each structure

    print '\m/ >.< \m/'
