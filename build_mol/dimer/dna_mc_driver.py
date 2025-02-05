#!/usr/bin/env python
#!/share/apps/bin/python
#
# Author:  Steven C. Howell
# Purpose: Run ds_dna_monte_carlo.py using gui-like input
# Created: 30 Octomber 2014
#
# $Id$
#
# 0000000011111111112222222222333333333344444444445555555555666666666677777777778
# 2345678901234567890123456789012345678901234567890123456789012345678901234567890


# import sassie.simulate.monte_carlo.ds_dna_monte_carlo.ds_dna_monte_carlo as ddmc
# import sassie.simulate.monte_carlo.ds_dna_monte_carlo.special.input_filter as input_filter
import dna.ds_dna_monte_carlo as ddmc
import dna.special.input_filter as input_filter

if __name__ == "__main__":

    svariables = {}

    # User Input
    svariables['runname'] = ('testing', 'string')
    svariables['path'] = ('./', 'string')
    svariables['infile'] = ('dimer.pdb', 'string')
    svariables['refpdb'] = ('dimer.pdb', 'string')
    svariables['psffile'] = ('dimer.psf', 'string')
    svariables['ofile'] = ('dimer.dcd', 'string')
    svariables['trials'] = ('100', 'int')
    svariables['goback'] = ('50', 'int')

    # Molecule Specific Input
    svariables['n_flex_regions'] = ('3', 'int')
    svariables['theta_max'] = ('10', 'float_array')
    # provide 's_theta_z_max' to use a different theta max for twisting vs
    # bending
    svariables['theta_z_max'] = ('1', 'float_array')
    svariables['first_res_per_region'] = ('6, , ', 'int_array')
    svariables['n_cont_res_per_region'] = ('14, , ', 'int_array')
    # s svariables['align_low_res']         = ('1', 'int')  # not yet implemented
    # s svariables['align_high_res']        = ('15', 'int')  # not yet
    # implemented
    svariables['dna_segnames'] = ('DNA1, DNA2', 'string_array')
    svariables['dna1_resids'] = ('5, end', 'int_array')
    svariables['dna2_resids'] = ('end-4, 5', 'int_array')
    svariables['n_rigid_groups'] = ('2', 'int')
    svariables['rigid_group1'] = (
        'A0, B0, C0, D0, E0, F0, G0, H0', 'string_array')
    svariables['rigid_group2'] = (
        'A1, B1, C1, D1, E1, F1, G1, H1', 'string_array')

    # Specialized/Advanced Inputs
    # set to N > 1 to have N base-pairs coarse-grained into 1 bead
    svariables['bp_per_bead'] = ('1', 'int')
    # set to N > 1 to apply rotations averaged over N coars-grained beads
    svariables['softrotation'] = ('1', 'int')
    # set to N > 1 to only record structures after every N trials
    svariables['n_dcd_write'] = ('1', 'int')
    svariables['seed'] = ('0', 'int')  # goback random seed
    # may be incorperated for electrostatics
    svariables['temperature'] = ('300', 'float')
    svariables['debug'] = ('False', 'bool')  # 'True' will display debug output
    # 'True' will generate a file labeling paird DNA base pairs for all flexible beads (useful for calculating dihedral angles)
    svariables['write_flex'] = ('True', 'bool')
    # 'True' will keep the coarse-grain DNA and protein pdb and dcd files
    svariables['keep_cg_files'] = ('False', 'bool')
    # 'False' will keep duplicate structure when move fails
    svariables['keep_unique'] = ('True', 'bool')
    # 'True' will remove the coarse-grain pkl file forcing a re-coarse-graining (can be time consuming often necessary)
    svariables['rm_pkl'] = ('True', 'bool')
    # 'True' will remove the coarse-grain pkl file forcing a re-coarse-graining (can be time consuming often necessary)
    svariables['openmm_min'] = ('False', 'bool')

    error, variables = input_filter.type_check_and_convert(svariables)
    vaiables = ddmc.prepare_dna_mc_input(variables)

    # Run
    ddmc.main(variables)
