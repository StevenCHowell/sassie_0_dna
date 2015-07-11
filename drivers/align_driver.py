import x_dna.drivers.myAlign as align
import os
'''
the inputs object should contain the following attributes
goal:    goal pdb
ref:     reference pdb containing molecule info for moving pdb/dcd
move:    pdb/dcd to align
out:     output dcd file
path:    output path
goal_filter:     goal basis filter
move_filter:     move basis filter

note: inputs.ref and inputs.move can often be the same pdb
'''

dimer_dcds = ['2x167_k010/run0/monte_carlo/2x167_k010_4d_1k.dcd',
            '2x167_k010/run1/monte_carlo/2x167_k010_10d_1k.dcd',
            '2x167_k010/run2/monte_carlo/2x167_k010_10d_5k.dcd',
            '2x167_k010/run3/monte_carlo/2x167_k010_10d_5k.dcd',
            '2x167_k010/run5/monte_carlo/2x167_k010_7d_5k.dcd',
            '2x167_k010/run4/monte_carlo/2x167_k010_7d_5k.dcd']

trimer_dcds = ['3x167_k010/run0/monte_carlo/3x167_k010_10d_1k.dcd',
               '3x167_k010/run1/monte_carlo/3x167_k010_4d_1k.dcd',
               '3x167_k010/run2/monte_carlo/3x167_k010_10d_5k.dcd',
               '3x167_k010/run3/monte_carlo/3x167_k010_10d_5k.dcd',
               '3x167_k010/run4/monte_carlo/3x167_k010_7d_5k.dcd',
               '3x167_k010/run5/monte_carlo/3x167_k010_7d_5k.dcd',
               '3x167_k050/run0/monte_carlo/3x167_k050_s1.dcd',
               '3x167_k050/run1/monte_carlo/3x167_k050_4d_1k.dcd',
               '3x167_k050/run2/monte_carlo/3x167_k050_10d_1k.dcd',
               '3x167_k050/run3/monte_carlo/3x167_k050_10d_5k.dcd',
               '3x167_k050/run4/monte_carlo/3x167_k050_10d_5k.dcd',
               '3x167_k050/run5/monte_carlo/3x167_k050_7d_5k.dcd',
               '3x167_k050/run6/monte_carlo/3x167_k050_7d_5k.dcd',
               '3x167_k100/run0/monte_carlo/3x167_k100_s1.dcd',
               '3x167_k100/run1/monte_carlo/3x167_k100_4d_1k.dcd',
               '3x167_k100/run2/monte_carlo/3x167_k100_10d_1k.dcd',
               '3x167_k100/run3/monte_carlo/3x167_k100_10d_5k.dcd',
               '3x167_k100/run4/monte_carlo/3x167_k100_10d_5k.dcd',
               '3x167_k100/run5/monte_carlo/3x167_k100_7d_5k.dcd',
               '3x167_k100/run6/monte_carlo/3x167_k100_7d_5k.dcd',
               '3x167_k100/run7/monte_carlo/3x167_k100_28bp_15d_10k.dcd',
               '3x167_k200/run0/monte_carlo/3x167_k200_4d_1k.dcd',
               '3x167_k200/run1/monte_carlo/3x167_k200_10d_1k.dcd',
               '3x167_k200/run2/monte_carlo/3x167_k200_10d_5k.dcd',
               '3x167_k200/run3/monte_carlo/3x167_k200_10d_5k.dcd',
               '3x167_k200/run4/monte_carlo/3x167_k200_7d_5k.dcd',
               '3x167_k200/run5/monte_carlo/3x167_k200_7d_5k.dcd',
               '3x167_k200/run6/monte_carlo/3x167_k200_6d_1k.dcd',
               '3x167_k200/run7/monte_carlo/3x167_k200_6d_1k.dcd',
               '3x167/run0/monte_carlo/3x167_10d_10k.dcd',
               '3x167/run10/monte_carlo/3x167_10d_10k.dcd',
               '3x167/run11/monte_carlo/3x167_10d_10k.dcd',
               '3x167/run1/monte_carlo/3x167_10d_10k.dcd',
               '3x167/run_1ncp/monte_carlo/3x167_28bp_10d_10k.dcd',
               '3x167/run28/monte_carlo/3x167_28bp_10d_10k.dcd',
               '3x167/run2/monte_carlo/3x167_10d_10k.dcd',
               '3x167/run3/monte_carlo/3x167_10d_10k.dcd',
               '3x167/run6_wco/monte_carlo/3x167_10d_1k.dcd',
               '3x167/run6_woco/monte_carlo/3x167_10d_1k.dcd',
               '3x167/run7_wco/monte_carlo/3x167_10d_1k.dcd',
               '3x167/run7_woco/monte_carlo/3x167_10d_1k.dcd',
               '3x167/run8/monte_carlo/3x167_10d_10k.dcd',
               '3x167/run9/monte_carlo/3x167_10d_10k.dcd']

tetramer_dcds = ['4x167_k010/run0/monte_carlo/4x167_k010_s1.dcd',
                 '4x167_k010/run1/monte_carlo/4x167_k010_10d_2k.dcd',
                 '4x167_k010/run2/monte_carlo/4x167_k010_4d_1k.dcd',
                 '4x167_k010/run3/monte_carlo/4x167_k010_10d_5k.dcd',
                 '4x167_k010/run4/monte_carlo/4x167_k010_10d_5k.dcd',
                 '4x167_k050/run0/monte_carlo/4x167_k050_4d_1k.dcd',
                 '4x167_k050/run1/monte_carlo/4x167_k050_10d_1k.dcd',
                 '4x167_k050/run2/monte_carlo/4x167_k050_10d_5k.dcd',
                 '4x167_k050/run3/monte_carlo/4x167_k050_10d_5k.dcd',
                 '4x167_k100/run0/monte_carlo/4x167_k100_s1.dcd',
                 '4x167_k100/run1/monte_carlo/4x167_k100_10d_8k.dcd',
                 '4x167_k100/run2/monte_carlo/4x167_k100_4d_1k.dcd',
                 '4x167_k100/run3/monte_carlo/4x167_k100_10d_5k.dcd',
                 '4x167_k100/run4/monte_carlo/4x167_k100_10d_5k.dcd',
                 '4x167_k100/run5/monte_carlo/4x167_k100_10d_2k.dcd',
                 '4x167_k100/run6/monte_carlo/4x167_k100_10d_1k.dcd',
                 '4x167_k100/run7/monte_carlo/4x167_k100_10d_1k.dcd',
                 '4x167_k100/run8/monte_carlo/4x167_k100_36bp_15d_20k.dcd',
                 '4x167_k100/run9/monte_carlo/4x167_k100_44bp_15d_20k.dcd',
                 '4x167_mg01/run0/monte_carlo/4x167_mg01_4d_1k.dcd',
                 '4x167_mg01/run1/monte_carlo/4x167_mg01_10d_1k.dcd',
                 '4x167_mg01/run2/monte_carlo/4x167_mg01_10d_5k.dcd',
                 '4x167_mg01/run3/monte_carlo/4x167_mg01_10d_5k.dcd',
                 '4x167_mg01/run4/monte_carlo/4x167_mg01_10d_2k.dcd',
                 '4x167_mg01/run5/monte_carlo/4x167_mg01_10d_1k.dcd',
                 '4x167_mg01/run6/monte_carlo/4x167_mg01_5d_20k.dcd',
                 '4x167_mg01/run7/monte_carlo/4x167_k010_10d_1k.dcd',
                 '4x167/run0/monte_carlo/4x167.dcd',
                 '4x167/run0/monte_carlo/test.dcd',
                 '4x167/run1/monte_carlo/4x167_5d_10k.dcd',
                 '4x167/run_1ncp/monte_carlo/4x167_28bp_10d_10k.dcd',
                 '4x167/run_28bp/monte_carlo/4x167_28bp_10d_10k.dcd',
                 '4x167/run2/monte_carlo/4x167_5d_10k.dcd',
                 '4x167/run3/monte_carlo/4x167_5d_10k.dcd',
                 '4x167/run4/monte_carlo/4x167_5d_20k.dcd',
                 '4x167/run5/monte_carlo/4x167_5d_20k.dcd',
                 '4x167/run6/monte_carlo/4x167_5d_20k.dcd',
                 '4x167/run7/monte_carlo/4x167_5d_20k.dcd']

start_dir = os.getcwd()

all_dcds = ['4x167/run0/monte_carlo/test.dcd']

dimer = '/home/schowell/myData/sassieRuns/2x167_k010_s1_min.pdb'
trimer = '/home/schowell/myData/sassieRuns/3x167_min.pdb'
tetramer = '/home/schowell/myData/sassieRuns/c11_min.pdb'



align_basis = ('((name[i] == "CA") and (segname[i] == "3H2A") and '
               '(resid[i] > 105) and (resid[i] < 115))')

inputs = align.inputs()
inputs.goal = tetramer
inputs.ref = tetramer

inputs.path = ''
inputs.goal_filter = align_basis
inputs.move_filter = align_basis

for dcd in all_dcds:
    inputs.move = os.path.join(start_dir, dcd)
    inputs.out = inputs.move.replace('.dcd', '_al.dcd')
    align.align(inputs)
    os.system('mv %s %s' % (inputs.out, inputs.move))


print '\m/ >.< \m/'