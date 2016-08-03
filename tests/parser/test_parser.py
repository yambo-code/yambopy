# Author:
# Tests for the yambopy library
# 
#
import unittest
import sys
import os
import argparse
import subprocess
import filecmp
from yamboparser import YamboFile, YamboFolder


class TestFolder(unittest.TestCase):
     
    def test_folder_list(self):
        fold = YamboFolder(os.path.dirname(os.path.realpath(__file__))+'/testdata')
        assert len (fold.yambofiles)==3


class TestFile(unittest.TestCase):

    def test_qp_parsing(self):
        fl = YamboFile('o-GW_run.10.720.qp',os.path.dirname(os.path.realpath(__file__))+'/testdata')
        assert len(fl.data.keys()) == 4  # more intelligent test needed
        assert  fl.type == 'output_gw'

    def test_l_parsing(self):
        fl = YamboFile('l-GW_run.8.480_em1d_ppa_HF_and_locXC_gw0_rim_cut_CPU_1',os.path.dirname(os.path.realpath(__file__))+'/testdata')
        assert  not fl.data 
        assert len(fl.warnings) ==1
        assert len(fl.errors) == 1
        assert  fl.type == 'log'

    def test_r_parsing(self):
        fl = YamboFile('r-GW_run.8.480_em1d_ppa_HF_and_locXC_gw0_rim_cut',os.path.dirname(os.path.realpath(__file__))+'/testdata')
        assert fl.type=='report'
        assert fl.kpoints
        assert not fl.data
