# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.bse.bse_absorption import YamboBSEAbsorptionSpectra
test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','bse')

class TestYamboBSEAbsorptionSpectra(unittest.TestCase):

    def test_yambobseabsorptionspeectra(self):

        #load databases
        lat  = YamboLatticeDB.from_db_file(os.path.join(test_path,'SAVE','ns.db1'))
        exc  = YamboExcitonDB.from_db_file(lat,folder=os.path.join(test_path,'yambo'))

        #open show analysis of bse absorption spectra
        bsea = YamboBSEAbsorptionSpectra(exc)
        bsea.plot()

        #bsea.write_json()

if __name__ == '__main__':
    unittest.main()
