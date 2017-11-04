from __future__ import print_function
#
# Copyright (C) 2017 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
import unittest
import os
import shutil as sh
from yambopy.analyse import YamboAnalyser
from yambopy.io.outputfile import YamboOut

test_path = os.path.join('..','..','tests','reference','gw_conv')

class TestYamboAnalyse(unittest.TestCase):
    def setUp(self):
        """ Read the yambo GW output files
        """
        print(test_path)
        sh.copytree(test_path,'gw_conv')
        for dirpath,dirnames,filenames in os.walk('gw_conv'):
            if YamboOut.has_output(dirpath):
                y = YamboOut(dirpath,save_folder='gw_conv')
                y.pack()

    def test_yamboanalyse_gw_si(self):
        """ Analyse the yambo GW .json output files
        """
        y = YamboAnalyser('gw_conv')
        print(y)
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
    
        #test get_path    
        #test path_plotting

    def tearDown(self):
        sh.rmtree('gw_conv') 
    
if __name__ == '__main__':
    unittest.main()
