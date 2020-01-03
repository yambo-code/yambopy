# Copyright (C) 2018 Henrique Pereira Coutada Miranda
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
from qepy.lattice import Path

test_path = os.path.join(os.path.dirname(__file__),'..','data','refs','gw_conv')

class TestYamboAnalyse(unittest.TestCase):
    def setUp(self):
        """ Read the yambo GW output files
        """
        if os.path.isdir('gw_conv'): sh.rmtree('gw_conv')
        sh.copytree(test_path,'gw_conv')
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

        #test getting data
        ks_bands,qp_bands = y.get_bands(tags='reference')
        ks_bands.plot(show=False)
        print(ks_bands)
        qp_bands.plot(show=False)
        print(qp_bands)

        #test get_path    
        path = Path([[[1.0,1.0,1.0],'G'],
                     [[0.0,0.5,0.5],'X'],
                     [[0.0,0.0,0.0],'G'],
                     [[0.5,0.0,0.0],'L']], [20,20,20])
        ks_bands_path,qp_bands_path = y.get_bands(tags='FFTGvecs',path_kpoints=path,type_calc=('gw'))
        qp_bands_path.plot(show=False) 

    def tearDown(self):
        sh.rmtree('gw_conv') 
    
if __name__ == '__main__':
    unittest.main()
