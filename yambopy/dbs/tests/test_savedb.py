# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
from yambopy.dbs.savedb import YamboSaveDB
from qepy.lattice import Path
test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','gw_conv')

class TestYamboSaveDB(unittest.TestCase):

    def test_yambosavedb(self):
        """ test savedb """

        #open savedb
        filename = os.path.join(test_path,'SAVE')
        ys = YamboSaveDB(filename)
        ys.get_fermi()
        ys.write_kpoints()
        str(ys)

    def tearDown(self):
        os.remove('kpts_full.dat')
        os.remove('kpts.dat')

if __name__ == '__main__':
    unittest.main()
