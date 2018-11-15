# Author: Henrique Pereira Coutada Miranda
# Tests for yambopy
# Si
#
from __future__ import print_function
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
this_file_path = os.path.dirname(__file__)
yambopy_path   = os.path.join(this_file_path,'..','..','scripts')
yambopy_script = os.path.join(yambopy_path,'yambopy')
ref_folder     = os.path.join(this_file_path,'..','data','refs','gw_conv')

class ItestYambopyGW(unittest.TestCase):
    def itest_yambopy_analysegw(self):
        """ Test the yambopy analysegw executable
        """
        os.system('{} analysegw {} FFTGvecs -bc 5 -kc 3 -bv 4 -kv 1 -nd'.format(yambopy_script,ref_folder))
        out = np.loadtxt('analyse_gw_conv/gw_conv_FFTGvecs.dat')
        ref = np.loadtxt('reference/si/analyse_gw_conv/gw_conv_FFTGvecs.dat')
        print("ref:")
        print(ref)
        print("out:")
        print(out)
        self.assertEqual(np.isclose(ref,out,atol=1e-3).all(),True)

        os.system('{} analysegw {} BndsRnXp -bc 5 -kc 3 -bv 4 -kv 1 -nd'.format(yambopy_script,ref_folder))
        out = np.loadtxt('analyse_gw_conv/gw_conv_BndsRnXp.dat')
        ref = np.loadtxt('reference/si/analyse_gw_conv/gw_conv_BndsRnXp.dat')
        print("ref:")
        print(ref)
        print("out:")
        print(out)
        self.assertEqual(np.isclose(ref,out,atol=1e-3).all(),True)


if __name__ == '__main__':
    # Count the number of errors
    nerrors = 0
    ul = unittest.TestLoader()
    tr = unittest.TextTestRunner(verbosity=2)

    #
    # Test GW on yambo
    #
    suite = ul.loadTestsFromTestCase(TestYambopyGW)
    nerrors += not tr.run(suite).wasSuccessful()
