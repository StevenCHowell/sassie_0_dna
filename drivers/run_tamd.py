#!/usr/bin/env python
#
# Auther: Steven C. Howell
# Purpose: Run TAMD in parrallel with consolidated output
# Created: 18 April 2016

import errno
import glob
import logging
import os
import shutil
import subprocess
import time

import multiprocessing as mp
import numpy as np
import os.path as op
import pandas as pd

import sasmol.sasmol as sasmol
# import sassie.util.basis_to_python as b2p

# LOGGER = logging.getLogger(__name__) #add module name manually

class cd:
    '''
    Context manager for changing the current working directory
    http://stackoverflow.com/questions/431684
    '''

    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class TamdInputs(object):

    def __init__(self, **kwargs):
        default_path = kwargs.get('default_path', './')
        self.charmm_exe = kwargs.get('charmm_exe',
                                     '/home/schowell/bin/charmm.exe')
        self.n_cpus = kwargs.get('n_cpus', 4)
        self.run_name = kwargs.get('run_name', 'tamd_run_0')
        self.dcd_fname = kwargs.get('dcd_fname', 'tamd_test.dcd')
        self.dcd_path = kwargs.get('dcd_path', default_path)
        self.inp_fname = kwargs.get('inp_fname', 'quick.inp')
        self.inp_path = kwargs.get('inp_path', default_path)
        self.pdb_fname = kwargs.get('pdb_fname', 'tamd_test.pdb')
        self.pdb_path = kwargs.get('pdb_path', default_path)
        logging.debug('charmm_exe: {}'.format(self.charmm_exe))
        logging.debug('n_cpus: {}'.format(self.n_cpus))
        logging.debug('run_name: {}'.format(self.run_name))
        logging.debug('dcd_fname: {}'.format(self.dcd_fname))
        logging.debug('dcd_path: {}'.format(self.dcd_path))
        logging.debug('inp_fname: {}'.format(self.inp_fname))
        logging.debug('inp_path: {}'.format(self.inp_path))
        logging.debug('pdb_fname: {}'.format(self.pdb_fname))
        logging.debug('pdb_path: {}'.format(self.pdb_path))

def append_bk(folder):
    if folder[-1] == '/':
        new_folder = op.split(folder)[0] + '_BK/'
    else:
        new_folder = folder + '_BK/'
    if op.exists(new_folder):
        append_bk(new_folder)
    shutil.move(folder, new_folder)
    print 'moved %s to %s' % (folder, new_folder)


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


def run_in_parallel(inputs):

    run_name = inputs.run_name
    dcd_fname = inputs.dcd_fname
    dcd_path = inputs.dcd_path
    pdb_fname = inputs.pdb_fname
    pdb_path = inputs.pdb_path
    inp_fname = inputs.inp_fname
    inp_path = inputs.inp_path

    tamd_exe = inputs.charmm_exe
    tamd_path = op.join(run_name, 'torsion_angle_md/')

    n_cpus = inputs.n_cpus

    if op.exists(tamd_path):
        logging.warning('moving existing run folder: %s'.format(tamd_path))
        append_bk(tamd_path)
        if tamd_path == dcd_path:
            dcd_path = op.split(dcd_path)[0] + '_BK/'
    else:
        print 'created foxs output folder: %s' % tamd_path
    mkdir_p(tamd_path)

    # check the input
    try:
        dcd_full_fname = op.join(dcd_path, dcd_fname)
        op.exists(dcd_full_fname)
    except:
        logging.error('No such file: {}'.format(dcd_full_fname))

    try:
        inp_full_fname = op.join(inp_path, inp_fname)
        op.exists(inp_full_fname)
    except:
        logging.error('No such file: {}'.format(inp_full_fname))

    try:
        pdb_full_fname = op.join(pdb_path, pdb_fname)
        op.exists(pdb_full_fname)
    except:
        logging.error('No such file: {}'.format(pdb_full_fname))

    try:
        op.exists(tamd_exe)
    except:
        logging.error('No such file: {}'.format(tamd_exe))

    # split the dcd into subfolders
    sub_dirs = split_dcd_to_pdbs(pdb_full_fname, dcd_full_fname, tamd_path)

    for sub_dir in sub_dirs:
        dst = op.join(sub_dir, inp_fname)
        shutil.copyfile(inp_full_fname, dst)

    n_runs = len(sub_dirs)
    n_runs_per_cpu = np.array([n_runs / n_cpus] * n_cpus, dtype=int)
    n_runs_per_cpu[:n_runs%n_cpus] += 1 # evenly spread the remainder

    # run tamd instance on each folder
    if n_cpus > 1:
        # setup the process list of lists
        process_list = []
        test_list  = []
        for i in xrange(max(n_runs_per_cpu)):
            test_list.append([])
            process_list.append([])

        for i_cpu in xrange(n_cpus):  # setup the processes
            for i_job in xrange(n_runs_per_cpu[i_cpu]):
                i_sub_dir = i_job * n_cpus + i_cpu
                test_list[i_job].append(i_sub_dir)
                tamd_args = (sub_dirs[i_sub_dir], inp_fname, tamd_exe)
                process_list[i_job].append(mp.Process(target=tamd, args=tamd_args))

        # run the process list of lists
        for processes in process_list:
            logging.info('starting {} TAMD runs'.format(len(processes)))
            for p in processes:  # start the processes
                p.start()

            for p in processes:  # exit the completed processes
                p.join()
            logging.info('finished {} TAMD runs'.format(len(processes)))

    else:
        for i in xrange(len(sub_dirs)):
            logging.info('starting TAMD run: {}'.format(i+1))
            tamd(sub_dirs[i], inp_fname, tamd_exe)

    # collect the output
    for sub_dir in sub_dirs:
        sub_files = glob.glob(op.join(sub_dir, '*'))
        split_dir = op.split(sub_dir)
        sub_index = split_dir[-1][-5:]
        dst_dir = split_dir[0]
        for sub_file in sub_files:
            fname = op.split(sub_file)[-1]
            if fname.lower()[-3:] == 'dcd':
                dst = op.join(dst_dir, 'tamd_dyn_{}.dcd'.format(sub_index))
            elif fname.lower()[-3:] == 'log':
                dst = op.join(dst_dir, 'min_{}.log'.format(sub_index))
            else:
                dst = op.join(dst_dir, fname)

            shutil.move(sub_file, dst)
        os.rmdir(sub_dir)

    shutil.copy(pdb_full_fname, dst_dir)

    logging.info('finished running TAMD')


def split_dcd_to_pdbs(pdb_full_fname, dcd_full_fname, starting_dir):

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_full_fname)

    dcd_file = mol.open_dcd_read(dcd_full_fname)
    total_frames = dcd_file[2]
    n_atoms = dcd_file[1]

    sub_dirs = []
    pdb_out_fname = 'temp_0.pdb'

    for i in xrange(total_frames):
        sub_dir = op.join(starting_dir, 'sub{}'.format(str(i+1).zfill(5)))
        sub_dirs.append(sub_dir)
        mkdir_p(sub_dir)

        with cd(sub_dir):
            mol.read_dcd_step(dcd_file, 0)
            mol.write_pdb(pdb_out_fname, 0, 'w')

    return sub_dirs


def tail(f, n=10):
    '''
    return the last n lines of f
    adapted from: http://stackoverflow.com/questions/136168
    '''
    tail_str = 'tail -n %s %s' % (str(n), f)
    stdin, stdout = os.popen2(tail_str)
    stdin.close()
    lines = stdout.readlines()
    stdout.close()
    return lines[:]


def tamd(sub_dir, inp_fname, tamd_exe):
    '''
    INPUT:  variable descriptions:
        sub_dir:    directory to process
        inp_fname:  the charmm run script
        tamd_exe:   the charmm executable used to run tamd

    OUTPUT:
        ...

    '''
    with cd(sub_dir):
        try:
            op.exists(inp_fname)
        except:
            logging.error('missing input file: {}'.format(op.join(os.getcwd(),
                                                                  inp_fname)))

        # run the calculation
        log_file = open('tamd.log', 'w')
        run_cmd = '{} -i {}'.format(tamd_exe, inp_fname)
        subprocess.call(run_cmd.split(), stderr=log_file, stdout=log_file)


if __name__ == '__main__':

    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

    # inputs = TamdInputs(default_path='script_test', inp_fname='quick.inp',
    # n_cpus=6)
    # inputs.run_name = op.join('script_test', inputs.run_name)

    # running with more cpus than dcd frames
    inputs = TamdInputs(inp_path='script_test', inp_fname='quick.inp',
                        n_cpus=12, dcd_fname='hiv1_gag_ab_cluster_30A.dcd',
                        pdb_path='script_test')

    # running with just one cpu
    # inputs = TamdInputs(inp_path='script_test', inp_fname='quick.inp',
                        # n_cpus=1, dcd_fname='hiv1_gag_ab_cluster_30A.dcd',
                        # pdb_path='script_test')

    run_in_parallel(inputs)

    logging.info('\m/ >.< \m/')

