#!/usr/bin/env python
#
# Author:  Steven C. Howell
# Purpose: to run crysol on the sassie generated structures
# Created: 18 March 2015
#
# $Id: $
#
#0000000011111111112222222222333333333344444444445555555555666666666677777777778
#2345678901234567890123456789012345678901234567890123456789012345678901234567890

import x_dna.drivers.my_crysol_driver as crysol

in_vars = crysol.inputs()
in_vars.runname = 'run1'
in_vars.dcdfile = 'dimer.dcd'
in_vars.dcdpath = 'run1/dna_mc/'
in_vars.pdbfile = 'dimer_mod.pdb'
in_vars.pdbpath = './'
# maxh ~= (Q_max x D_max)/pi
in_vars.maxh = '18' # Dimer_maxh ~= (0.2 x 280)/pi
# in_vars.maxh = '21' # Trimer_maxh ~= (0.2 x 325)/pi     
# in_vars.maxh = '21' # Tetramer_maxh ~= (0.2 x 325)/pi     
in_vars.fib  = '18' # more does not slow crysol down noticeably
in_vars.maxs = '0.2'
in_vars.numpoints = '100'

o = crysol.Drv()
o.run_me(in_vars)