# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
from yambopy.bse.bse_absorption import YamboBSEAbsorptionSpectra
test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','bse')

class TestYamboBSEAbsorptionSpectra(unittest.TestCase):

    #@unittest.skip("Creating flows to handle tasks and use ypp inside a flow folder")
    def test_yambobseabsorptionspeectra(self):

        pass
        #open qpdb
        #bsea = YamboBSEAbsorptionSpectra('bse',path=test_path)
        #print(bsea)

        #bsea.get_excitons()

if __name__ == '__main__':
    unittest.main()
