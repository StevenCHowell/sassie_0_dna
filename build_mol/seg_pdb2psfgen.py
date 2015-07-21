#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Replace DNA sequence with another sequence
# Created: 24 April 2014
#
# $Id$
#
'''
This script loads a pdb structure file of DNA, and creates a '*.patches' file
with the psfgen patches needed to use psfgen to create the structure.
After running this script, the patches can be pasted into a psfgen file.
'''

import sassie.sasmol.sasmol as sasmol
import sys
import string
import time
import logging


def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        # prog='',
        # usage='',
        description='test functionality of the cgDNA move module',
        #epilog = 'no epilog found'
    )

    parser.add_argument("-p", "--pdb", help="all atom pdb file")
    parser.add_argument("-s", "--segnames", nargs='+',
                        help="segnames to extract")

    return parser.parse_args()


def main(inputs):
    mol_all = sasmol.SasMol(0)
    mol_all.read_pdb(inputs.pdb)

    print inputs.segnames

    segname1 = inputs.segnames[0]
    segname2 = inputs.segnames[1]

    print 'segname 1: ', segname1
    print 'segname 2: ', segname2

    psfgenFile = inputs.pdb[:-4] + '_patches.txt'

    outfile = open(psfgenFile, 'w')      # open the file
    timestr = time.strftime("# created on %d %B %Y by 'pdb2psfgen.py'\n")
    outfile.write(timestr)
    outfile.write('# dna1: segname ' + segname1 + '\n')
    outfile.write('# dna2: segname ' + segname2 + '\n')

    mol_segs = sasmol.SasMol(0)
    basis_filter = "( (segname[i]=='%s') or (segname[i]=='%s') )" % (
        segname1, segname2)
    error, mask = mol_all.get_subset_mask(basis_filter)
    error = mol_all.copy_molecule_using_mask(mol_segs, mask, 0)

    names = mol_segs.resname()
    ids = mol_segs.resid()
    c = mol_segs.segname()

    pyr = ['C', 'T', 'DC', 'DT', 'CYT', 'THY']
    pur = ['A', 'G', 'DA', 'DG', 'ADE', 'GUA']
    pyrStr = 'patch DEO1 '
    purStr = 'patch DEO2 '

    n = 0
    for (j, i) in enumerate(ids):
        # only want this to happend once for each residue
        if n != i:
            n = i
            # print 'adding line %d' % i
            skip = False
            if c[j] in segname1:
                dna = 'dna1:%d\n' % i
            elif c[j] in segname2:
                dna = 'dna2:%d\n' % i
            else:
                print 'Skipping residue from unspecified segname: ', c[j]
                skip = True

            if not skip:
                if names[j] in pyr:
                    outfile.write(pyrStr + dna)
                    # print pyrStr + dna
                elif names[j] in pur:
                    outfile.write(purStr + dna)
                    # print purStr + dna
                else:
                    print 'ERROR!!! unknown resname in specified segname: ', names[j]
                    print '\n'

    outfile.close()

    print 'COMPLETE \m/ >.< \m/'

if __name__ == "__main__":

    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' % __name__, level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    # make ARGS global
    ARGS = parse()
    main(ARGS)
