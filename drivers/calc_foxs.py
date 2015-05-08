#!/usr/bin/env python
# Auther: Steven C. Howell
# Purpose: Run crysol in parrallel with consolidated output
# Created: 10/09/2014
# $Id: parallel_crysol.py 183 2015-03-19 22:42:34Z xraylab $

import sassie.sasmol.sasmol as sasmol
import numpy as np
import os.path as op
import glob
import multiprocessing as mp

import sys
import os
import errno
import subprocess
import logging
import cmd
import shutil
import time


def foxs(res_dir, dcd_name, first_last, pdb_name, foxs_exe, output, max_q=0.2, 
         num_points=50):
    '''
    FOXS is the function to read in structures from a DCD/PDB file and 
    calculate SAXS profiles using the binary program foxs.exe

    INPUT:  variable descriptions:

        res_dir:    directory to move results to
        pdb_name:   name of pdb file for loading atom info
        dcd_name:   name of dcd file for loading atom coordinates
        foxs_exe:   foxs executable path+filename
        output:     output object
        max_q:      pmaximum q value (default = 0.2)
        num_points: number of points in the profile

    OUTPUT:

        output.

        files stored in res_dir/ directory:
        run_name*.dat:        SAXS profile from crysol (q -vs I(q))

    REFERENCE:
    D. Schneidman-Duhovny, et al. Biophysical Journal 2013. 105 (4), 962-974
    D. Schneidman-Duhovny, et al. NAR 2010. 38 Suppl:W540-4 
    '''

    
    # start multiple runs to calculate scattering in each subfolder
    # each run should move the output back to the original foxs folder

    tempdir=run_name+'/crysol/'+run_name+'_files'
    tempdir='./tmp/'

    m1=sasmol.SasMol(0)
    m1.read_pdb(pdb_full_name)
    check_atom_names(m1)

    try:
        if(dcdfile[-3:] == 'dcd'):
            ldcdfile = m1.open_dcd_read(dcdpath+'/'+dcdfile)
            nf = ldcdfile[2]
            intype = 'dcd'
        elif(dcdfile[-3:] =='pdb'):
            m1.read_pdb(dcdpath+'/'+dcdfile)
            check_atom_names(m1)
            nf = m1.number_of_frames()
            intype = 'pdb'
    except:
        message='input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
        message+=' :  stopping here'
        print_failure(message,txtOutput)

    print 'number of frames = ',nf
    print 'intype = ',intype   

    # first create file name arrays

    pdb=[] ; newpdb=[]
    for i in range(nf):
        nst = str(i+1).zfill(5) 	
        pdb.append(pfile1+'_'+nst)

    print 'there are %d files to process\n' % (nf)

    k=0

    # get OS type (assumes cryson27.osx is for MAC & crysol25.l86 for linux

    osst = 'uname -a'
    value = os.popen(osst,'r').readlines()
    svalue = string.split(value[0])
    os_type = svalue[0]

    if(os_type != "Darwin"):
        os_type = "Linux"

    print 'os_type = ',os_type

    ttxt=time.ctime()
    st=''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" %(st))
    txtOutput.put("DATA FROM RUN: %s \n\n" %(ttxt))

    # now loop over files and process them using crysol
    numproc=0	
    for i in range(nf):
        numproc+=1
        k=k+1
        n=i+1;
        currfilebase=pdb[i]
        currfile=tempdir+'/dum.pdb'
        if(intype == 'dcd'):
            m1.read_dcd_step(ldcdfile,i)
            m1.write_pdb(currfile,0,'w') ; time.sleep(0.25)
        elif(intype == 'pdb'):
            m1.write_pdb(currfile,i,'w') ; time.sleep(0.25)
        if(os_type == "Linux"):
            currfile=tempdir+'/dum.pdb'
        print currfile+'\n'
        currcrys=tempdir+'/'+currfilebase+'.ans'
        curroutfile=open(currcrys,mode='w+')
        if(os_type == "Linux"):
            curroutfile.write("%s\n" % (currfile))
            curroutfile.write("%s\n" % (option))
            curroutfile.write(currfile+"\n")
        elif(os_type == "Darwin"):
            curroutfile.write("%s\n" % (option))
            curroutfile.write("%s\n" % (currfile))
        curroutfile.write("%s\n" % (maxh))
        curroutfile.write("%s\n" % (fib))
        curroutfile.write("%s\n" % (maxs))
        curroutfile.write("%s\n" % (numpoints))
        curroutfile.write("%s\n" % (hydrogens))
        curroutfile.write("N\n")
        curroutfile.write("%s\n" % (edensolv))
        curroutfile.write("%s\n" % (contrast))
        curroutfile.write("\n") # accept default
        curroutfile.write("\n") # accept default
        curroutfile.close()
        currcommand=cryexe+'<' + currcrys
        os.system(currcommand)
        if(delafs==0):
            os.system('mv -f *.alm ' + foxs_path+'/'+currfilebase+'.alm')
            os.system('mv -f *.flm ' + foxs_path+'/'+currfilebase+'.flm')
            os.system('mv -f *.sav ' + foxs_path+'/'+currfilebase+'.sav')   
        else:
            os.system('rm -f *.alm')
            os.system('rm -f *.flm')
            os.system('rm -f *.sav')
        if k==1:
            os.system('mv '+tempdir+'/*.ans ' + foxs_path+'/')
        else:
            os.system('rm -f ' +tempdir+'/*.ans')
        mvst='mv  *.log ' +foxs_path+'/'+currfilebase+'.log'
        print mvst
        os.system(mvst)
        mvst='mv  *.int ' +foxs_path+'/'+currfilebase+'.int'
        print mvst
        os.system(mvst)
        os.system('rm -f '+currfile)

        fraction_done = (float(i+1)/float(nf))
        progress_string='COMPLETED '+str(i+1)+' of '+str(nf)+' : '+str(fraction_done*100.0)+' % done'
        print('%s\n' % progress_string)
        print('%s\n' % progress_string)
        report_string='STATUS\t'+str(fraction_done)
        txtOutput.put(report_string)

    os.system('rm -Rf '+tempdir)

    if(intype == 'dcd'):
        m1.close_dcd_read(ldcdfile[0])

    txtOutput.put("Processed %s DCD frames\n" % (numproc))
    txtOutput.put("Data stored in directory: %s\n\n" % ('./'+foxs_path))
    txtOutput.put("\n%s \n" %(st))
    time.sleep(0.5)

    output.put('FoXS finished calculating %d profiles in %s' % (n)



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
    else:
        shutil.move(folder, new_folder)
        print 'moved %s to %s' % (folder, new_folder)

class folder_exists(cmd.Cmd):

    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = '(0/1/2)> '

    def do_move(self, arg):
        append_bk(self.runname)
        return True

    def help_move(self):
        print '-- move run folder to run_BK'

    def do_replace(self, arg):
        print 'removing run folder'
        shutil.rmtree(self.runname)
        return True

    def help_replace(self):
        print 'remove and replace run folder'

    def do_quit(self, arg):
        print 'exiting program'
        sys.exit(1)

    def help_quit(self):
        print '-- terminates the application'

    def default(self, arg):
        print 'invalid selection, please select: 0/1/2'


    #shortcuts
    do_0 = do_quit
    do_1 = do_move
    do_2 = do_replace
    help_0 = help_quit
    help_1 = help_move
    help_2 = help_replace

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

def collect_crysol(sub_dirs, runname, sleep):
    out_dir = os.getcwd() + '/' + runname + '/crysol'
    mkdir_p(out_dir)

    n_out_files = 1

    for (i, sub_dir) in enumerate(sub_dirs):
        logging.debug('waiting for %s' % sub_dir)
        with cd(sub_dir):
            # read the end of the file test.out
            # out_file = 'par_crysol_%02d.out' % (i+1)
            out_files = glob.glob('*out')
            not_done = True
            while not_done:
                time.sleep(sleep)
                for (j, out_file) in enumerate(out_files):
                    if tail(out_file, 1) == ['GCRYSOL IS DONE\n']:
                        not_done = False

            print '%s is complete' % sub_dir
            sub_files = glob.glob('*/crysol/*int')
            more_sub_files = glob.glob('crysol/*int')
            n_sub_files = len(sub_files)
            n_more_sub_files = len(more_sub_files)
            if n_sub_files > 0 and n_more_sub_files > 0:
                logging.warning(('found crysol output in both "./*/crysol/" '
                                 '(%d files) and "./crysol/" (%d files), order may be '
                                 'unexpected') % (n_sub_files, n_more_sub_files) )
            for another_file in more_sub_files:
                sub_files.append(another_file)

            error = sub_files.sort()
            for (j, sub_file) in enumerate(sub_files):
                logging.debug('moving %s' % sub_file)
                file_name = sub_file[:-4]
                new_name = runname + '_' + str(n_out_files + j).zfill(5)
                logging.debug('moving %s%s.int to %s/%s.int' % (sub_dir, 
                                                                file_name, out_dir, new_name))
                os.system('mv %s.int %s/%s.int' % (file_name, out_dir, new_name))
                os.system('mv %s.log %s/%s.log' % (file_name, out_dir, new_name))                

            # os.system('mv *out *dcd *pdb ../')
            os.system('mv *out ../')
            n_out_files += len(sub_files)

def split_dcd(pdb_full_name, dcd_full_name, n_cpus, starting_dir):

    mol = sasmol.SasMol(0)
    mol.read_pdb(pdb_full_name)

    dcd_file = mol.open_dcd_read(dcd_full_name)
    total_frames = dcd_file[2]
    n_atoms = dcd_file[1]
    # copy_mask = np.ones(n_atoms, dtype=np.int32)
    copy_mask = mol.get_subset_mask('all')

    n_frames_sub = total_frames/n_cpus
    last_frame = 0
    sub_dirs = []
    sub_dcd_names = []
    first_last = []
    for cpu in xrange(1, n_cpus+1):
        sub_dir = op.join(starting_dir, 'sub%s' % str(cpu).zfill(2))
        sub_dirs.append(sub_dir)
        mkdir_p(sub_dir)
        os.system('cp %s %s' % pdb_full_name, sub_dir)
        sub_mol = sasmol.SasMol(0)
        mol.copy_molecule_using_mask(sub_mol, copy_mask, 0)
        with cd(sub_dir):
            if cpu == n_cpus:
                n_frames_sub = n_frames_sub + total_frames % n_cpus
            dcd_out_name = 'sub%s.dcd' % str(cpu).zfill(2)
            sub_dcd_names.append(dcd_out_name)
            first = last_frame
            last = last_frame + n_frames_sub
            dcd_out_file = sub_mol.open_dcd_write(dcd_out_name)
            for (i, frame) in enumerate(xrange(first, last)):
                sub_mol.read_dcd_step(dcd_file, frame)
                sub_mol.write_dcd_step(dcd_out_file, 0, i+1)

            sub_mol.close_dcd_write(dcd_out_file)

        first_last.append([first, last])
        last_frame += n_frames_sub

    return sub_dirs, sub_dcd_names, first_last

def main(inputs):

    run_name   = inputs.run_name 
    dcd_name   = inputs.dcd_name 
    pdb_name   = inputs.pdb_name 
    pdb_path   = inputs.pdb_path 
    dcd_path   = inputs.dcd_path 
    
    foxs_exe   = inputs.foxs_exe  
    max_q      = inputs.max_q    
    num_points = inputs.num_points

    n_cpus     = inputs.n_cpus

    pdb_full_name = inputs.pdb_full_name = op.join(pdb_path, pdb_name)
    dcd_full_name = inputs.dcd_full_name = op.join(dcd_path, dcd_name)

    #check the input
    assert op.exists(dcd_full_name), 'ERROR: no such file "%s"' % dcd_full_name 
    assert op.exists(pdb_full_name), 'ERROR: no such file "%s"' % pdb_full_name 
    assert op.exists(foxs_exe), 'ERROR: no such file "%s"' % foxs_exe

    foxs_path = inputs.foxs_path = op.join(run_name, 'foxs')
    if op.exists(fox_path):
        print 'WARNING: run folder exists (%s), moving it\n' % fox_path
        append_bk(fox_path)
    else:
        print 'created foxs output folder: %s' % fox_path
    mkdir_p(foxs_path)

    # split the dcd into subfolders
    sub_dirs, sub_dcd_names, first_last = split_dcd(
        pdb_full_name, dcd_full_name, n_cpus, foxs_path)
    
    # run foxs instance on each folder
    output = mp.Queue()
    processes = []
    for cpu in n_cpus:  # setup the processes
        foxs_args = (sub_dirs[cpu], sub_dcd_names[cpu], first_last[cpu], 
                     pdb_name, foxs_exe, max_q, num_points, output)
        processes.append(mp.Process(target=foxs, args=foxs_args))

    for p in processes: # start the processes
        p.start()

    for p in processes: # exit the completed processes
        p.join()

    results = [output.get() for p in processes]
    #collect the results
    collect_crysol(sub_dirs, inputs.runname, inputs.sleep)

    print 'finished %d foxs calculations\n     \m/ >.< \m/ ' % first_last[-1][1]

if __name__ == '__main__':
    import argparse

    if '-v' in sys.argv:
        logging.basicConfig(filename='%s.log' %__file__[:-3], level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    ARGS = parse()    
    main(ARGS)
