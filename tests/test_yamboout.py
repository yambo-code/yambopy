from __future__ import print_function
#
# Author: Henrique Pereira Coutada Miranda
# Tests for yambopy
# Si
#
import matplotlib
import unittest
import sys
import os
import shutil
import argparse
import subprocess
import filecmp
import shutil as sh
from yambopy import *
from qepy import *

class TestYamboOut(unittest.TestCase):
    """ This class creates the input files for Si and compares them to reference files
    """
    def test_yamboout(self):

        yo = YamboOut('reference/gw')
        print(yo)
        assert yo.logs == ['l-yambo_em1d_ppa_HF_and_locXC_gw0_CPU_1', 
                           'l-yambo_em1d_ppa_HF_and_locXC_gw0_CPU_2']
        assert yo.run == ['r-yambo_em1d_ppa_HF_and_locXC_gw0'] 
        assert yo.output == ['o-yambo.qp'] 
        assert yo.netcdf == ['ndb.HF_and_locXC','ndb.QP'] 

if __name__ == '__main__':
    # Count the number of errors
    nerrors = 0
    ul = unittest.TestLoader()
    tr = unittest.TextTestRunner(verbosity=2)

    #
    # Test pw.x
    #
    suite = ul.loadTestsFromTestCase(TestYamboOut)
    nerrors += not tr.run(suite).wasSuccessful()

    sys.exit(nerrors)
