#!/usr/bin/env python
# Author:  Steven C. Howell
# Purpose: Prepare PDB for modeling
# Created: 04/24/2014
# $Id: makeHistones.py,v 1.3 2014-08-07 21:02:14 schowell Exp $
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

class residue:
    def __init__(self, resid, resname):
        self.resid = resid
        self.resname = resname
    def __repr__(self):
        return repr((self.resid, self.resname))
        

def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        #prog='',
        #usage='',
        description = 'test functionality of the cgDNA move module',
        #epilog = 'no epilog found'
    )

    parser.add_argument("-p", "--pdb", help="all atom pdb file")
    parser.add_argument("-c", "--chains", help="chains to extract, e.g.: [\"\'x\'\",\"\'y\'\"] or [\\'x\\',\\'y\\']")

    return parser.parse_args()

def main():
    # aa_pdb = '../1zbb_tetra_uncombined.pdb'
    # aa_pdb = '1zbb_original.pdb'

    aa = sasmol.SasMol(0)
    aa.read_pdb(ARGS.pdb)
    chain_mols = []
    errors = []


    amino_acids = {'ALA': 'A',
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
                   'VAL': 'V'}
    
    dna = {'DG': 'G',
           'DA': 'A',
           'DT': 'T',
           'DC': 'C', 
           'GUA': 'G', 
           'ADE': 'A', 
           'THY': 'T', 
           'CYT': 'C'}
    
    chains = ['I','J']
    print 'chains =', chains
    print 'ARGS.chains =', ARGS.chains
    
    
    for chain in chains:
        if chain.lower() == chain:
            chain_name = 'chain_' + chain + '0'
        else:
            chain_name = 'chain_' + chain + '1'
        chain_name = ARGS.pdb[:-4] + chain_name
        basis_filter = "(chain[i] == '" + chain + "')"
        error, mask = aa.get_subset_mask(basis_filter)
        if error: print error
        chain_mol = sasmol.SasMol(0)
        error = aa.copy_molecule_using_mask(chain_mol, mask, 0)
        if error: print error
        chain_mol.write_pdb(chain_name+'.pdb', 0, 'w')
        chain_mols.append(chain_mol)
    
        # resids.sort()
        resA = 0
        res_min = np.min(chain_mol.resids())
        res_max = np.max(chain_mol.resids())
        print 'min resid:', res_min
        
        resA = 0
        residue_list = []
        
        ## create a sorted list of the residues
        for (i, resB) in enumerate(chain_mol.resid()):
            if resB != resA:
                # print 'chain_mol.resname()[i]:', chain_mol.resname()[i]
                # print 'chain_mol.resid()[i]:', chain_mol.resid()[i]
                residue_list.append(residue(resB, chain_mol.resname()[i]))
                resA = resB
            
        residue_sequence = sorted(residue_list, key=lambda residue: residue.resid)
        #with open(chain_name+'.txt', 'w') as outFile:
            #for res in residue_sequence:
                #outFile.write(str(res.resid) + '\t' + res.resname + '\n')
    
        with open(chain_name+'.seq', 'w') as outFile:          
            if chain_mol.moltypes() == ['protein']:
                for (i, res) in enumerate(residue_sequence):
                    outFile.write(amino_acids[res.resname])
                    if 0 == (i + 1) % 50:
                        outFile.write('\n')
            elif chain_mol.moltypes() == ['dna']:
                for (i, res) in enumerate(residue_sequence):
                    outFile.write(dna[res.resname])
                    if 0 == (i + 1) % 50:
                        outFile.write('\n')
            else:
                print 'ERROR, unexpected molecule type'
            
        #s for resB in resids:
        #s     if resB != resA + 1:
        #s         print 'missing residue/s, skipped chain', chain, 'btwn residues:', resA, resB
        #s     resA = resB
        print 'max resid:', res_max
        
        print 'finished chain', chain_name
    
    
    print 'COMPLETE :(|) '

if __name__ == "__main__":

    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' %__name__, level=logging.DEBUG)
    else:
        logging.basicConfig()

    # make ARGS global
    ARGS = parse()
    main()

    
    
