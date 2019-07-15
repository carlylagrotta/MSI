import sys
sys.path.append('..')

from cti_core import soln2cti as cti_pp

import cantera as ct
gas = ct.Solution("../data/test_data/lam.cti")

cti_pp.write(gas, "../data/test_data/cti_test.cti")
