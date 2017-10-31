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

        #pack the data
        yo.pack()

        #verify ndb.QP data
        keys = sorted(yo.files['ndb.QP'].keys())
        assert keys == ['Band', 'E', 'E-Eo', 'Eo', 'Kpoint', 'Kpoint_index', 'Z', 'qp_table', 'type']

        #verify o-*.qp data
        keys = sorted(yo.files['o-yambo.qp'].keys())
        assert keys == ['Band', 'E-Eo', 'Eo', 'K-point', 'Sc|Eo', 'input', 'type']

        #compare data from ndb.QP and o-yambo.qp
        of = yo.files['o-yambo.qp']['Eo']
        qo = yo.files['ndb.QP']['Eo']
        #the difference is larger than 1e-3
        assert np.all(np.isclose(of,qo,atol=5e-3))

class TestYamboAnalyse(unittest.TestCase):
    def test_yamboout_gw_si(self):
        """ Read the yambo GW output files
        """
        sh.rmtree('gw_conv')
        sh.copytree('reference/gw_conv','gw_conv')
        for dirpath,dirnames,filenames in os.walk('gw_conv'):
            if YamboOut.has_output(dirpath):
                y = YamboOut(dirpath,save_folder='gw_conv')
                y.pack()

    def test_yamboanalyse_gw_si(self):
        """ Analyse the yambo GW .json output files
        """
        y = YamboAnalyser('gw_conv')
        netcdf_files = y.get_files_type('netcdf_gw')
        keys = sorted(netcdf_files.keys())
        assert keys == ['BndsRnXp_1_20', 'BndsRnXp_1_30', 'FFTGvecs_00010', 'FFTGvecs_00015', 
                        'NGsBlkXp_00002', 'NGsBlkXp_00005', 'reference']

        netcdf_files = y.get_files_type('netcdf_gw','FFTGvecs')
        keys = sorted(netcdf_files.keys())
        assert keys == ['FFTGvecs_00010', 'FFTGvecs_00015']

        netcdf_files = y.get_files_type('netcdf_gw',('FFTGvecs','reference'))
        keys = sorted(netcdf_files.keys())
        assert keys == ['FFTGvecs_00010', 'FFTGvecs_00015','reference']

        #test plotting
        gw_bands = y.plot_gw(tags='FFTGvecs')
        print(gw_bands)
        ks_bands = y.plot_ks(tags='reference')
        print(ks_bands)

if __name__ == '__main__':
    # Count the number of errors
    nerrors = 0
    ul = unittest.TestLoader()
    tr = unittest.TextTestRunner(verbosity=2)

    suite = ul.loadTestsFromTestCase(TestYamboOut)
    nerrors += not tr.run(suite).wasSuccessful()

    suite = ul.loadTestsFromTestCase(TestYamboAnalyse)
    nerrors += not tr.run(suite).wasSuccessful()


    sys.exit(nerrors)
