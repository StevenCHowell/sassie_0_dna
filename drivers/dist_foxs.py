#!/usr/bin/env python
# Auther: Steven C. Howell
# Purpose: Check that there is no overlap then run crysol (parrallel)
# Created: 10/09/2014
# $Id: parallel_crysol.py 183 2015-03-19 22:42:34Z xraylab $

import glob, shutil, os, subprocess, time
import sassie.sasmol.sasmol as sasmol
import numpy as np
import os.path as op
import multiprocessing as mp
import x_dna.util.basis_to_python as b2p
import x_dna.ds_dna_monte_carlo as ddmc
import pandas as pd
# import sys
# import errno
# import logging
# import cmd

class inputs():
    def __init__(self, parent = None):
        pass

def foxs(sub_dir, dcd_name, first_last, pdb_full_name, foxs_exe, basis='all', 
         max_q=0.2, num_points=50):
    '''
    FOXS is the function to read in structures from a DCD/PDB file and 
    calculate SAXS profiles using the binary program foxs.exe

    INPUT:  variable descriptions:

        sub_dir:    directory to move results to
        dcd_name:   name of dcd file for loading atom coordinates
        first_last: index of the dcd frames used in the calculation
        pdb_full_name:   name of pdb file for loading atom info
        foxs_exe:   foxs executable path+filename
        output:     output object
        basis:      basis for calculating the scattering (default = all) 
        max_q:      pmaximum q value (default = 0.2)
        num_points: number of points in the profile

    OUTPUT:

        files stored in res_dir/ directory:
        run_name*.dat:        SAXS profile from crysol (q -vs I(q))

    REFERENCE:
    D. Schneidman-Duhovny, et al. Biophysical Journal 2013. 105 (4), 962-974
    D. Schneidman-Duhovny, et al. NAR 2010. 38 Suppl:W540-4 
    '''

    
    # start multiple runs to calculate scattering in each subfolder
    # each run should move the output back to the original foxs folder

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_full_name)

    # setup the molecule for calculating scattering using FoXS
    _, mask = mol.get_subset_mask(basis)
    assert sum(mask) > 0, 'ERROR: no atom selected using basis: %s' % basis
    sub_mol = sasmol.SasMol(0)
    mol.copy_molecule_using_mask(sub_mol, mask, 0)

    # setup the molecule for checking for DNA overlap
    p_basis = ' name[i] == "P" '
    _, p_mask = mol.get_subset_mask(p_basis)
    assert sum(p_mask) > 0, 'ERROR: no atoms selected using basis: %s' % p_basis
    p_mol = sasmol.SasMol(0)
    mol.copy_molecule_using_mask(p_mol, p_mask, 0)

    # setup the log file
    path, name = op.split(sub_dir)
    name += '.log'
    log_file = open(op.join(path,name), 'w')
    with cd(sub_dir):
        assert op.exists(dcd_name), 'ERROR: invalid input file:%s' % dcd_name
        dcd_file = mol.open_dcd_read(dcd_name)
        nf = dcd_file[2]
        log_file.write('Beginning the FoXS calculation for %d files\n' % (nf))
        ten_percent = nf/10 # this is intentionally an int
        
        # create file name arrays
        res_dir = op.split(sub_dir)[0]
        file_prefix = op.split(res_dir)[1]
        out_base = [file_prefix + '_' + str(i+1).zfill(5) for i in 
                    xrange(first_last[0], first_last[1])]
    
        # setup the DCD and Rg output
        dcd_out_fname = '../%s_filtered.dcd' % file_prefix
        dcd_out_file = mol.open_dcd_write(dcd_out_fname)
        n_written = 0
        rg = []
        label = []
        overlapped = []
        
        st = ''.join(['=' for x in xrange(60)])
        log_file.write(("\n%s \n" %(st)))
        time_txt = time.ctime()
        log_file.write("DATA FROM RUN: %s \n\n" %(time_txt))
        # now loop over files and process them using foxs
        tic = time.time()
        
        for i_proc in xrange(nf):
            n_proc = i_proc + 1
            
            # setup the input file
            mol.read_dcd_step(dcd_file, i_proc)
            p_mol.setCoor(mol.get_coor_using_mask(0, p_mask)[1])
            cutoff = 3
            if not ddmc.f_overlap1(p_mol.coor()[0], cutoff):
                n_written += 1
                mol.write_dcd_step(dcd_out_file, 0, n_written)
                rg.append(mol.calcrg(0))
                label.append(out_base[n_written-1]) # index from 1 (not 0)
        
                cur_file = op.join(out_base[n_written] + '.pdb')
                sub_mol.setCoor(mol.get_coor_using_mask(0, mask)[1])
                sub_mol.write_pdb(cur_file, 0, 'w')
    
                # run the calculation
                foxs_out = subprocess.check_output(
                    ['foxs', '-q', str(max_q), '-s', str(num_points), cur_file], 
                    stderr=subprocess.STDOUT)
                log_file.write(foxs_out)
                
                # move and cleanup the output
                os.remove(cur_file)              # remove the pdb file
                os.remove(cur_file[:-3] + 'plt') # remove the plt file
                tmp_out = cur_file + '.dat'
                cur_out = '../' + cur_file[:-3] + 'dat'
                os.rename(tmp_out, cur_out)      # move the I(q) file
            
                if 0 == np.mod(n_proc, ten_percent):
                    log_file.write('moved %s to %s\n' % (tmp_out, cur_out))
                    fraction_done = n_proc/float(nf) * 100
                    log_file.write('COMPLETED %s, %d of %d: %0.0f %% done\n' % (
                        cur_file, n_proc, nf, fraction_done))
            else:
                overlapped.append(n_proc)
        toc = time.time() - tic
        rg_dict = {'rg': rg}
        rg_df = pd.DataFrame(rg_dict, index=label)
        rg_df.index.name = 'labels'
        rg_df.to_csv('../rg.csv', sep='\t')
        mol.close_dcd_write(dcd_out_file)
        mol.close_dcd_read(dcd_file[0])
        os.remove(dcd_name)
        try:
            os.remove(atomid.out)
        except:
            pass
    os.rmdir(sub_dir)
    
    log_file.write("Data stored in directory: %s\n\n" % res_dir)
    log_file.write("FoXS calculated %s DCD frames in %d s (%d frames/min)\n" 
                   % (n_written, toc, int(n_written/toc*60)))
    print "FoXS calculated %s DCD frames in %d s (%d frames/min)\n" % (
        n_written, toc, int(n_written/toc*60))
    log_file.close()

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

def append_bk(folder):
    new_folder = folder + '_BK'
    if op.exists(new_folder):
        append_bk(new_folder)    
    shutil.move(folder, new_folder)
    print 'moved %s to %s' % (folder, new_folder)

class cd:
    """
    Context manager for changing the current working directory
    http://stackoverflow.com/questions/431684
    """
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def tail(f, n=10):
    '''
    return the last n lines of f
    adapted from: http://stackoverflow.com/questions/136168
    '''
    tail_str = 'tail -n %s %s' % (str(n), f)
    stdin,stdout = os.popen2(tail_str)
    stdin.close()
    lines = stdout.readlines()
    stdout.close()
    return lines[:]

def split_dcd(pdb_full_name, dcd_full_name, n_cpus, starting_dir):

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_full_name)

    dcd_file = mol.open_dcd_read(dcd_full_name)
    total_frames = dcd_file[2]
    n_atoms = dcd_file[1]
    # copy_mask = np.ones(n_atoms, dtype=np.int32)
    _, copy_mask = mol.get_subset_mask('all')

    n_frames_sub = total_frames/n_cpus
    last_frame = 0
    sub_dirs = []
    sub_dcd_names = []
    first_last = []
    for cpu in xrange(1, n_cpus+1):
        sub_dir = op.join(starting_dir, 'sub%s' % str(cpu).zfill(2))
        sub_dirs.append(sub_dir)
        mkdir_p(sub_dir)
        sub_mol = sasmol.SasMol(0)
        mol.copy_molecule_using_mask(sub_mol, copy_mask, 0)
        with cd(sub_dir):
            if cpu == n_cpus:
                n_frames_sub = n_frames_sub + total_frames % n_cpus
            dcd_out_name = 'sub%s.dcd' % str(cpu).zfill(2)
            sub_dcd_names.append(dcd_out_name)
            first = last_frame
            last = last_frame + n_frames_sub
            if n_cpus == 1:
                rel_dcd_name = '../../../%s'%dcd_full_name
                assert op.exists(rel_dcd_name), 'ERROR: did not find dcd file'
                subprocess.call(['ln', '-s', rel_dcd_name, dcd_out_name])
            else:
                dcd_out_file = sub_mol.open_dcd_write(dcd_out_name)
                for (i, frame) in enumerate(xrange(first, last)):
                    sub_mol.read_dcd_step(dcd_file, frame)
                    sub_mol.write_dcd_step(dcd_out_file, 0, i+1)

                sub_mol.close_dcd_write(dcd_out_file)

        first_last.append([first, last])
        last_frame += n_frames_sub
    mol.close_dcd_read(dcd_file[0])
    
    return sub_dirs, sub_dcd_names, first_last

def main(inputs):

    run_name   = inputs.run_name 
    dcd_name   = inputs.dcd_name 
    dcd_path   = inputs.dcd_path 
    pdb_name   = inputs.pdb_name 
    pdb_path   = inputs.pdb_path 
    
    foxs_exe   = inputs.foxs_exe  
    try:
        max_q = inputs.max_q    
    except:    
        max_q = 0.2
    try:
        num_points = inputs.num_points
    except:
        num_points = 50
    try:
        basis = inputs.basis
    except:
        basis = 'all'
    # n_cpus     = inputs.n_cpus # not fully implemented for parallel
    n_cpus = 1
    
    pdb_full_name = inputs.pdb_full_name = op.join(pdb_path, pdb_name)
    dcd_full_name = inputs.dcd_full_name = op.join(dcd_path, dcd_name)

    #check the input
    assert op.exists(dcd_full_name), 'ERROR: no such file "%s"' % dcd_full_name 
    assert op.exists(pdb_full_name), 'ERROR: no such file "%s"' % pdb_full_name 
    assert op.exists(foxs_exe), 'ERROR: no such file "%s"' % foxs_exe

    foxs_path = inputs.foxs_path = op.join(run_name, 'foxs')
    if op.exists(foxs_path):
        print 'WARNING: run folder exists (%s), moving it\n' % foxs_path
        append_bk(foxs_path)
    else:
        print 'created foxs output folder: %s' % foxs_path
    mkdir_p(foxs_path)

    # run foxs instance on each folder
    py_basis = b2p.parse_basis(basis)
    # split the dcd into subfolders
    sub_dirs, sub_dcd_names, first_last = split_dcd(
        pdb_full_name, dcd_full_name, n_cpus, foxs_path)        
    processes = []
    
    if n_cpus > 1:
        for cpu in xrange(n_cpus):  # setup the processes
            foxs_args = (sub_dirs[cpu], sub_dcd_names[cpu], first_last[cpu], 
                         pdb_full_name, foxs_exe, py_basis, max_q, num_points)
            processes.append(mp.Process(target=foxs, args=foxs_args) )
            
    
        for p in processes: # start the processes
            p.start()
    
        for p in processes: # exit the completed processes
            p.join()
    else:
        foxs(sub_dirs[0], sub_dcd_names[0], first_last[0], pdb_full_name, 
             foxs_exe, py_basis, max_q, num_points)

    print 'finished %d foxs calculations\n' % first_last[-1][1]

if __name__ == '__main__':
    # main(ARGS)
    NotImplemented
