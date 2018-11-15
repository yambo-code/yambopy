# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function
import unittest
import os
import shutil as sh
import numpy as np
from yambopy.io.inputfile import YamboIn

test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','gw_conv')

class TestYamboIn(unittest.TestCase):
    """ This class creates the input files for Si and compares them to reference files
    """
    def test_yamboin(self):

        #read input file
        yi = YamboIn.from_file(folder=test_path)
        ref_file = str(yi)

        #write input file
        os.makedirs('optimize')
        yi.write('optimize/yambo.in')

        #read again
        yi = YamboIn.from_file(folder='optimize')
        assert len(str(yi)) == len(ref_file)

        #test convergence
        conv = { 'FFTGvecs': [[5,10,15],'Ry'],
                 'NGsBlkXp': [[1,2,5], 'Ry'],
                 'BndsRnXp': [[1,10],[1,20],[1,30]] }
        yi.folder = 'optimize'
        yi.optimize(conv,folder='optimize')

    def tearDown(self):
        sh.rmtree('optimize')

if __name__ == "__main__":
    unittest.main() 
