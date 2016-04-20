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
import sassie.sasmol.sasmol as sasmol
import numpy as np


def main():
    ''' preprocessing for commandline executions'''
    import logging
    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' % __name__, level=logging.DEBUG)
    else:
        logging.basicConfig()

    parser = argparse.ArgumentParser(  # prog='',  #usage='',
        description='Separate the PDB into individual PDB and Sequence files by segname',
    )

    parser.add_argument("pdbfile", help="file of pdb structure")
    parser.add_argument(
        "-o", "--outfile", nargs='?', help="file to save sequence")
    args = parser.parse_args()

    return pdb_get_sequence(args.pdbfile, args.outfile)


def pdb_get_sequence(pdbobj=None, outfile=None):
    ''' get the sequence of a sasmol object '''

    if isinstance(pdbobj, basestring):
        pdbfile = pdbobj
        pdbobj = sasmol.SasMol(0)
        pdbobj.read_pdb(pdbfile)

    resname2seq = {'ALA': 'A',  # amino acids
                   'ARG': 'R',
                   'ASN': 'N',
                   'ASP': 'D',
                   'CYS': 'C',
                   'GLU': 'E',
                   'GLN': 'Q',
                   'GLY': 'G',
                   'HIS': 'H',
                   'ILE': 'I',
                   'LEU': 'L',
                   'LYS': 'K',
                   'MET': 'M',
                   'PHE': 'F',
                   'PRO': 'P',
                   'SER': 'S',
                   'THR': 'T',
                   'TRP': 'W',
                   'TYR': 'Y',
                   'VAL': 'V',
                   'HSE': 'H',
                   'G': 'G',    # DNA
                   'A': 'A',
                   'T': 'T',
                   'C': 'C',
                   'DG': 'G',
                   'DA': 'A',
                   'DT': 'T',
                   'DC': 'C',
                   'GUA': 'G',
                   'ADE': 'A',
                   'THY': 'T',
                   'CYT': 'C'}

    resid_all = pdbobj.resid()
    idx_unique = np.nonzero(np.insert(resid_all[1:] - resid_all[0:-1], 0, 1))
    idx_unique = idx_unique[0]  # it appears to be a tuple

    resname_all = pdbobj.resname()
    sequence = map(lambda i: resname2seq[resname_all[i]], idx_unique)

    if outfile == None:
        print "Sequence: total {} residues".format(len(sequence))
        print "".join(sequence)
    else:
        with open(outfile, 'w') as fileobj:
            fileobj.write("".join(sequence))

    return sequence

if __name__ == "__main__":
    main()
