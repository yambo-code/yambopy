# Author: 
# Tests for the yambopy library
# 
#
from __future__ import print_function
import unittest
import sys
import os
import argparse
import subprocess
import filecmp
from yamboparser import YamboFile, YamboFolder

folder = os.path.join(os.path.dirname(__file__),'..','..','yambopy','data','refs','parser')

class TestFolder(unittest.TestCase):
     
    def test_folder_list(self):
        fold = YamboFolder(os.path.join(folder,'t2_parse_qps/'))
        assert len (fold.yambofiles)==5


class TestFileT1(unittest.TestCase):

    def test_qp_parsing(self):
        fl = YamboFile('o-GW_run.10.720.qp',os.path.join(folder,'t1_errors_warnings'))
        assert len(list(fl.data.keys())) == 4  # more intelligent test needed
        assert  fl.type == 'output_gw'

    def test_l_parsing(self):
        fl = YamboFile('l-GW_run.8.480_em1d_ppa_HF_and_locXC_gw0_rim_cut_CPU_1',os.path.join(folder,'t1_errors_warnings'))
        assert  not fl.data 
        assert len(fl.warnings) ==1
        assert len(fl.errors) == 1
        assert  fl.type == 'log'

    def test_r_parsing(self):
        fl = YamboFile('r-GW_run.8.480_em1d_ppa_HF_and_locXC_gw0_rim_cut',os.path.join(folder,'t1_errors_warnings'))
        assert fl.type=='report'
        assert fl.kpoints
        assert not fl.data

class TestFileT2(unittest.TestCase):

    def test_qp_parsing(self):
        fl = YamboFile('o-yambo.qp',os.path.join(folder,'t2_parse_qps'))
        assert  fl.type == 'output_gw'

    def test_l_parsing(self):
        fl = YamboFile('l-yambo_em1d_HF_and_locXC_gw0',os.path.join(folder,'t2_parse_qps'))
        assert  fl.type == 'log'

    def test_r_parsing(self):
        fl = YamboFile('r-yambo_em1d_life',os.path.join(folder,'t2_parse_qps'))
        assert fl.type=='report'

        fl = YamboFile('r-yambo_em1d_HF_and_locXC_gw0',os.path.join(folder,'t2_parse_qps'))
        assert fl.type=='report'

    def test_ndb_qp_parsing(self):
        fl = YamboFile('ndb.QP',os.path.join(folder,'t3_parse_netcdf'))
        print("fl type", fl.type)
        assert fl.type=='netcdf_gw'

    def test_ndb_hf_parsing(self):
        fl = YamboFile('ndb.HF_and_locXC',os.path.join(folder,'t3_parse_netcdf'))
        print("fl type", fl.type)
        assert fl.type=='netcdf_hf'

if __name__ == "__main__":

    #t1_errors_warnings
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFileT1)
    unittest.TextTestRunner(verbosity=2).run(suite)

    #t2_parse_qps
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFileT2)
    unittest.TextTestRunner(verbosity=2).run(suite)

