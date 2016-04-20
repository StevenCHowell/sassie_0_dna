'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import sys
import string
import locale
import bisect
import time
import numpy
import sassie.sasmol.sasmol as sasmol

#       ALIGN2
#
#	12/17/2004	--	adapted from trehalose project for gag modeling :	jc
#	12/20/2004	--	align N-terminal of CA in 1L6N and 1E6J		:	jc
#	10/16/2005	--	generic align structures 			:	jc
#	03/11/2009	--	modified to allow for different sized molecules	:	jc
#	01/01/2010	--	added sasmol support				:	jc
#	07/29/2011	--	adapted to use sasop & sasmath classes		:	jc
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789

'''
        ALIGN2 is the module that overlaps molecules from a dcd/pdb file
	onto another molecule over a given basis.  The two molecule types
	do not need to be the same but the number of basis atoms used for
	the overlap do need to be identical.

        This module is called from Align Frames from the main
        GUI through the graphical_align.py script.

	REFERENCE:

	W. Kabsch
    	Acta Crystallog. sect. A  32  922-923  (1976)

    	W. Kabsch
    	Acta Crystallog. sect. A  34  827-828  (1978)

'''


def print_failure(message, txtOutput):

    txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n")
    txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
    txtOutput.put(message)

    return


def unpack_variables(variables):

    runname = variables['runname'][0]
    path = variables['path'][0]
    infile = variables['infile'][0]
    pdbmol1 = variables['pdbmol1'][0]
    pdbmol2 = variables['pdbmol2'][0]
    ofile = variables['ofile'][0]
    basis1 = variables['basis1'][0]
    basis2 = variables['basis2'][0]
    lowres1 = variables['lowres1'][0]
    lowres2 = variables['lowres2'][0]
    highres1 = variables['highres1'][0]
    highres2 = variables['highres2'][0]

    return runname, path, infile, pdbmol1, pdbmol2, basis1, lowres1, highres1, basis2, lowres2, highres2, ofile


def align(variables, txtOutput):
    '''
    ALIGN is the function to read in variables from GUI input and
    overlap the molecules in a dcd/pdb file onto the coordinates of
    a reference pdb structure over a given basis.

            runname: 	        project name
            path:                   input/output filepath
            pdbmol1:                reference pdb (mol 1)
            pdbmol2:                input pdb file (mol 2)
            infile:                 input (pdb or dcd) filename (mol 2)
            basis1:                 basis for molecule 1
            basis2:                 basis for molecule 2
            lowres1:                low residue for overlap molecule 1
            highres1:               high residue for overlap molecule 1
            lowres2:                low residue for overlap molecule 2
            highres2:               high residue for overlap molecule 2

    OUTPUT:

            files stored in "runname"/align directory

            ofile:			output filename
            ofile*.minmax:		text file with min & max dimensions

    '''

    runname, path, infile, pdbmol1, pdbmol2, basis1, lowres1, highres1, basis2, lowres2, highres2, ofile = unpack_variables(
        variables)

    alignpath = runname + '/align/'
    direxist = os.path.exists(alignpath)
    if(direxist == 0):
        os.system('mkdir -p ' + alignpath)

    print 'runname = ', runname

    dcd = []
    dcd.append(infile)
    ndcd = 1
    minmaxfile = ofile + '.minmax'
    mmfile = open(alignpath + minmaxfile, 'w')

    ttxt = time.ctime()

    st = ''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" % (st))
    txtOutput.put("DATA FROM RUN: %s \n\n" % (ttxt))

    m1 = sasmol.SasMol(0)
    m2 = sasmol.SasMol(1)

    m1.readpdb(path + pdbmol1)
    m2.readpdb(path + pdbmol2)

    try:
        if(infile[-3:] == 'dcd'):
            m2.readdcd(path + infile)

        elif(infile[-3:] == 'pdb'):
            m2.readpdb(path + infile)

    except:
        message = 'input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
        message += ' :  stopping here'
        print_failure(message, txtOutput)

    nf2 = m2.number_of_frames()

    txtOutput.put("Total number of frames = %d\n\n" % (nf2))

    mass1 = m1.mass()
    mass2 = m2.mass()

    name1 = m1.name()
    name2 = m2.name()

    basis_filter_1 = 'name[i] == "' + basis1 + '" and (resid[i] >= ' + str(
        lowres1) + ' and resid[i] <= ' + str(highres1) + ')'
    basis_filter_2 = 'name[i] == "' + basis2 + '" and (resid[i] >= ' + str(
        lowres2) + ' and resid[i] <= ' + str(highres2) + ')'

    error, mask1 = m1.get_subset_mask(basis_filter_1)
    error, mask2 = m2.get_subset_mask(basis_filter_2)

    print 'numpy.sum(mask1) = ', numpy.sum(mask1)
    print 'numpy.sum(mask2) = ', numpy.sum(mask2)

    sub_m1 = sasmol.SasMol(2)
    error = m1.copy_molecule_using_mask(sub_m1, mask1, 0)
    print 'error = ', error

    sub_m2 = sasmol.SasMol(3)
    error = m2.copy_molecule_using_mask(sub_m2, mask2, 0)
    print 'error = ', error

    com_sub_m1 = sub_m1.calccom(0)
    sub_m1.center(0)
    coor_sub_m1 = sub_m1.coor()[0]

    print 'com_sub_m1 = ', com_sub_m1

    for i in xrange(nf2):

        m2.center(i)

        error, sub_m2.coor = m2.get_coor_using_mask(i, mask2)
        sub_m2.setCoor(sub_m2.coor)
        com_sub_m2 = sub_m2.calccom(0)
        sub_m2.center(0)
        coor_sub_m2 = sub_m2.coor[0]

        m2.align(i, coor_sub_m2, com_sub_m2, coor_sub_m1, com_sub_m1)

        if(((i + 1) % (float(nf2) / 10.0) == 0 or (nf2 < 10))):
            fraction_done = (float(i + 1) / float(nf2))
            progress_string = 'COMPLETED ' + \
                str(i + 1) + ' of ' + str(nf2) + ' : ' + \
                str(fraction_done * 100.0) + ' % done'
            print('%s\n' % progress_string)
            report_string = 'STATUS\t' + str(fraction_done)
            txtOutput.put(report_string)

    try:
        if(ofile[-3:] == 'dcd'):
            print ' writing DCD file'
            m2.writedcd(alignpath + ofile)

        elif(ofile[-3:] == 'pdb' and nf2 == 1):
            print ' writing PDB file'
            m2.writepdb(alignpath + ofile, 0, 'w')
        elif(ofile[-3:] == 'pdb' and nf2 > 1):
            print ' writing PDB file'
            for i in xrange(nf2):
                if(i == 0):
                    m2.writepdb(alignpath + ofile, i, 'w')
                else:
                    m2.writepdb(alignpath + ofile, i, 'a')

        else:
            message = 'output filename ' + ofile + \
                ' needs to end in either ".pdb" (1 frame) or ".dcd" (1 or more frames)\n'
            message += ' :  writing output file as a ' + ofile + '.dcd\n'
            print '\n\n', message, '\n\n'
            print ' writing DCD file'
            ofile = ofile + '.dcd'
            m2.writedcd(alignpath + ofile)

    except:
        message = 'Could not write output file'
        print_failure(message, txtOutput)

    total_min_array, total_max_array = m2.calcminmax()

    min_x = total_min_array[0]
    max_x = total_max_array[0]
    min_y = total_min_array[1]
    max_y = total_max_array[1]
    min_z = total_min_array[2]
    max_z = total_max_array[2]

    txtOutput.put("minimum x = %lf\t maximum x = %lf -> range: %lf Angstroms\n" %
                  (min_x, max_x, (max_x - min_x)))
    txtOutput.put("minimum y = %lf\t maximum y = %lf -> range: %lf Angstroms\n" %
                  (min_y, max_y, (max_y - min_y)))
    txtOutput.put("minimum z = %lf\t maximum z = %lf -> range: %lf Angstroms\n\n" %
                  (min_z, max_z, (max_z - min_z)))

    print 'Aligned data (nf=%i) were written to %s\n' % (nf2, './' + alignpath + ofile)
    txtOutput.put("\nAligned data (nf=%i) were written to %s\n\n" %
                  (nf2, './' + alignpath + ofile))
    txtOutput.put("\n%s \n" % (st))
    time.sleep(0.5)

    print 'ALIGN2 IS DONE'
    return()
