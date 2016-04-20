#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Prepare PDB for modeling
# Created: 24 April 2014
#

'''
This script creates a separate pdb for each chain of 'aa_pdb'
It also creates a sequence file for the residue sequence of that chain
'''

import sys
import os
import os.path as op
import subprocess
import logging
import sassie.sasmol.sasmol as sasmol
import numpy as np
import pdb_get_sequence


def main():
    ''' preprocessing for commandline executions'''
    import logging
    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' % __name__, level=logging.DEBUG)
    else:
        logging.basicConfig()

    parser = argparse.ArgumentParser(  # prog='', #usage='',
        description='Separate the PDB into individual PDB and Sequence files by segname',
    )

    parser.add_argument("-p", "--pdbfile", help="all atom pdb file")
    parser.add_argument(
        "-o", "--outfile", nargs='?', help="prefix of outfile to save")
    parser.add_argument("-g", "--get_seq", default=True,
                        action='store_true', help="whether get sequence")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-s", "--segnames", nargs='+', help="segnames to extract")
    group.add_argument(
        "-c", "--chainids", nargs='+', help="segnames to extract")

    args = parser.parse_args()
    return pdb_get_chains(args.pdbfile, outfile=args.outfile,
                          segnames=args.segnames, chainids=args.chainids, get_seq=args.get_seq)


def pdb_get_chains(pdbobj=None, outfile='seg_', segnames=None, chainids=None, get_seq=True):
    ''' get the sequence of a sasmol object '''

    # filename is passed, read it
    if isinstance(pdbobj, basestring):
        pdbfile = pdbobj
        pdbobj = sasmol.SasMol(0)
        pdbobj.read_pdb(pdbfile)

    # set the filter prefix (chainid has priority)
    if segnames != None:
        filter_name = segnames
        filter_tmpl = "(segname[i] == '{}')"
    if chainids != None:
        filter_name = chainids
        filter_tmpl = "(chain[i] == '{}')"

    # a single string is passed, convert it to a list
    if isinstance(filter_name, basestring):
        filter_name = [filter_name]

    seg_mols = []
    for eachfilter in filter_name:
        print "Filter pdb by: ", filter_tmpl.format(eachfilter)
        error, mask = pdbobj.get_subset_mask(filter_tmpl.format(eachfilter))
        if error:
            print error

        eachfilter_mol = sasmol.SasMol(0)
        error = pdbobj.copy_molecule_using_mask(eachfilter_mol, mask, 0)
        if error:
            print error

        # eachfilter_mol.setSegname(eachfilter)
        eachfilter_mol.write_pdb(outfile + '.pdb', 0, 'w')
        seg_mols.append(eachfilter_mol)

        if get_seq:
            pdb_get_sequence.pdb_get_sequence(eachfilter_mol)

        print 'COMPLETE'
    return seg_mols

if __name__ == "__main__":
    main()
