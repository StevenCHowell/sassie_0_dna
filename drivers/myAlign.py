#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: Align several structures using an align basis
# Created: 6 February 2015

'''
This driver serves as a wrapper tool the align method used to align a
pdb/dcd to the coordinates of a goal pdb structure.
'''

import logging
import sys
import os
import sassie.sasmol.sasmol as sasmol
# import sassie.tools.align2 as a2
# import sassie_0_dna.util.align2 as a2



class inputs():

    def __init__(self, parent=None):
        pass


def parse():
    ''' Returns arguments in parser'''

    parser = argparse.ArgumentParser(
        # prog='',
        # usage='',
        description='test functionality of the cgDNA move module',
        #epilog = 'no epilog found'
    )

    parser.add_argument("-g", "--goal", help="goal pdb")
    parser.add_argument("-r", "--ref",
                        help="pdb with atom info for pdb/dcd that is moving")
    parser.add_argument("-m", "--move", help="pdb/dcd to be moved/aligned")
    parser.add_argument("-o", "--out", help="output dcd file name")
    parser.add_argument("-ms", "--move_seg_chain",
                        help="segname or chain to match")
    parser.add_argument("-gs", "--goal_seg_chain",
                        help="segname or chain to match")
    parser.add_argument("-msc", "--move_seg_or_chain",
                        help="matching a segname or chain")
    parser.add_argument("-gsc", "--goal_seg_or_chain",
                        help="matching a segname or chain")
    parser.add_argument("-p", "--path", help="output path")
    parser.add_argument("-gmn", "--goal_min",
                        help="minimun residue to match on goal molecule")
    parser.add_argument("-gmx", "--goal_max",
                        help="maximum residue to match on goal molecule")
    parser.add_argument("-mmn", "--move_min",
                        help="minimun residue to match on move molecule")
    parser.add_argument("-mmx", "--move_max",
                        help="miximum residue to match on move molecule")
    parser.add_argument("-ba", "--basis_atoms", help="basis_atoms to match")

    return parser.parse_args()


def align(inputs):
    '''
    input:
    ------
        inputs: object should contain the following attributes
            goal:    goal pdb
            ref:     reference pdb containing molecule info for moving pdb/dcd
            move:    pdb/dcd to align
            out:     output dcd file
            path:    output path
            goal_filter:     goal basis filter
            move_filter:     move basis filter

    note: inputs.ref and inputs.move can ofter be the same pdb
    '''
    aa_goal_pdb = inputs.goal
    aa_move_pdb = inputs.ref
    aa_move_file = inputs.move
    save_file = inputs.out
    path = inputs.path

    try:
        goal_filter = inputs.goal_filter
    except:
        basis_atoms = inputs.basis_atoms
        goal_seg_or_ch = inputs.goal_seg_or_chain
        goal_segname = inputs.goal_seg_chain
        goal_res_max = inputs.goal_max
        goal_res_min = inputs.goal_min
        goal_filter = ('((%s[i] == "%s") and (name[i] == "%s") and '
                       '(resid[i] >= %s) and (resid[i] <= %s))' % (
                           goal_seg_or_ch, goal_segname, basis_atoms,
                           goal_res_min, goal_res_max))

    try:
        move_filter = inputs.move_filter
    except:
        basis_atoms = inputs.basis_atoms
        move_seg_or_ch = inputs.move_seg_or_chain
        move_segname = inputs.move_seg_chain
        move_res_max = inputs.move_max
        move_res_min = inputs.move_min
        move_filter = ('((%s[i] == "%s") and (name[i] == "%s") and '
                       '(resid[i] >= %s) and (resid[i] <= %s))' % (
                           move_seg_or_ch, move_segname, basis_atoms,
                           move_res_min, move_res_max))

    # check input
    assert os.path.exists(aa_move_file), ('ERROR: no such file - %s' %
                                          aa_move_file)
    assert os.path.exists(aa_move_pdb), ('ERROR: no such file - %s' %
                                         aa_move_pdb)
    assert os.path.exists(aa_goal_pdb), ('ERROR: no such file - %s' %
                                         aa_goal_pdb)

    # create the SasMol objects
    sub_goal = sasmol.SasMol(0)
    sub_move = sasmol.SasMol(0)
    aa_goal = sasmol.SasMol(0)
    aa_move = sasmol.SasMol(0)

    aa_goal.read_pdb(aa_goal_pdb)
    aa_move.read_pdb(aa_move_pdb)

    if aa_move_file[-3:] == 'pdb':
        aa_move.read_pdb(aa_move_file)
        n_frames = aa_move.number_of_frames()
        in_type = 'pdb'
    elif aa_move_file[-3:] == 'dcd':
        dcd_file = aa_move.open_dcd_read(aa_move_file)
        n_frames = dcd_file[2]
        in_type = 'dcd'
    else:
        message = "\n~~~ ERROR, unknown input type ~~~\n"
        print_failure(message, txtOutput)
        return

    out_type = save_file[-3:].lower()
    if 'dcd' == out_type:
        dcd_out_file = aa_move.open_dcd_write(path + save_file)
    elif 'pdb' == out_type:
        dcd_out_file = None

    error, goal_seg_mask = aa_goal.get_subset_mask(goal_filter)
    assert not error, error
    error, move_seg_mask = aa_move.get_subset_mask(move_filter)
    assert not error, error

    error = aa_goal.copy_molecule_using_mask(sub_goal, goal_seg_mask, 0)
    assert not error, error
    error = aa_move.copy_molecule_using_mask(sub_move, move_seg_mask, 0)
    assert not error, error

    # calculate the center of mass of the subset of m1
    com_sub_goal = sub_goal.calccom(0)
    sub_goal.center(0)                         # center the m1 coordinates
    # get the m1 centered coordinates
    coor_sub_goal = sub_goal.coor()[0]

    for i in xrange(n_frames):
        if in_type == 'dcd':
            aa_move.read_dcd_step(dcd_file, i)
            # move m2 to be centered at the origin
            aa_move.center(0)
            error, sub_move.coor = aa_move.get_coor_using_mask(
                0, move_seg_mask)
            sub_move.setCoor(sub_move.coor)
            # calculate the center of mass of the subset of m2
            com_sub_move = sub_move.calccom(0)
            # move the subset of m2 to be centered at the origin
            sub_move.center(0)
            # get the new coordinates of the subset of m2
            coor_sub_move = sub_move.coor[0]
            # align m2 using the transformation from sub_m2 to sub_m1
            aa_move.align(
                0, coor_sub_move, com_sub_move, coor_sub_goal, com_sub_goal)
        elif in_type == 'pdb':
            # move m2 to be centered at the origin
            aa_move.center(i)
            error, sub_move.coor = aa_move.get_coor_using_mask(
                i, move_seg_mask)
            sub_move.setCoor(sub_move.coor)
            # calculate the center of mass of the subset of m2
            com_sub_move = sub_move.calccom(0)
            # move the subset of m2 to be centered at the origin
            sub_move.center(0)
            # get the new coordinates of the subset of m2
            coor_sub_move = sub_move.coor[0]
            # align m2 using the transformation from sub_m2 to sub_m1
            aa_move.align(
                i, coor_sub_move, com_sub_move, coor_sub_goal, com_sub_goal)

        aa_move.write_dcd_step(dcd_out_file, 0, i + 1)

    if in_type == 'dcd':
        aa_move.close_dcd_read(dcd_file[0])
    if out_type == 'dcd':
        aa_move.close_dcd_write(dcd_out_file)

    print 'COMPLETE \m/ >.< \m/'


def align_mol(inputs):
    '''
    input:
    ------
        intputs: object that should contain the following attributes
            aa_goal:    goal sasmol object
            aa_move:    sasmol object to align
            goal_basis: goal basis for alignment
            move_basis: move basis for alignment

    returns:
    --------
        out: aligned sasmol object

    note: inputs.ref and inputs.move are typically the same pdb/dcd
    '''
    aa_goal = inputs.aa_goal
    aa_move = inputs.aa_move
    goal_basis = inputs.goal_basis
    move_basis = inputs.move_basis

    # create the SasMol objects
    sub_goal = sasmol.SasMol(0)
    sub_move = sasmol.SasMol(0)

    error, goal_seg_mask = aa_goal.get_subset_mask(goal_basis)
    error, move_seg_mask = aa_move.get_subset_mask(move_basis)

    error = aa_goal.copy_molecule_using_mask(sub_goal, goal_seg_mask, 0)
    error = aa_move.copy_molecule_using_mask(sub_move, move_seg_mask, 0)

    # calculate the center of mass of the subset of m1
    com_sub_goal = sub_goal.calccom(0)
    sub_goal.center(0)                         # center the m1 coordinates
    # get the m1 centered coordinates
    coor_sub_goal = sub_goal.coor()[0]

    aa_move.center(0)                  # move m2 to be centered at the origin
    error, sub_move.coor = aa_move.get_coor_using_mask(0, move_seg_mask)
    sub_move.setCoor(sub_move.coor)
    # calculate the center of mass of the subset of m2
    com_sub_move = sub_move.calccom(0)
    # move the subset of m2 to be centered at the origin
    sub_move.center(0)
    # get the new coordinates of the subset of m2
    coor_sub_move = sub_move.coor[0]
    # align m2 using the transformation from sub_m2 to sub_m1
    aa_move.align(0, coor_sub_move, com_sub_move, coor_sub_goal, com_sub_goal)

    # the return statement may not be necessary
    # return aa_move

if __name__ == "__main__":

    import argparse
    if '-v' in sys.argv:
        logging.basicConfig(filename='_log-%s' % __name__, level=logging.DEBUG)
        sys.argv.pop(sys.argv.index('-v'))
    else:
        logging.basicConfig()

    # make ARGS global
    ARGS = parse()

    align(ARGS)
