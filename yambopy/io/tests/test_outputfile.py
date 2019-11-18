# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function
import unittest
import os
import numpy as np
from yambopy.io.outputfile import YamboOut

test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','gw')

class TestYamboOut(unittest.TestCase):
    """ This class creates the input files for Si and compares them to reference files
    """
    def test_yamboout(self):

        yo = YamboOut(test_path)
        print(yo)
        assert yo.logs == ['l-yambo_em1d_ppa_HF_and_locXC_gw0_CPU_1',
                           'l-yambo_em1d_ppa_HF_and_locXC_gw0_CPU_2']
        assert yo.run == ['r-yambo_em1d_ppa_HF_and_locXC_gw0']
        assert yo.output == ['o-yambo.qp']
        assert sorted(yo.netcdf) == sorted(['ndb.HF_and_locXC','ndb.QP'])

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

if __name__ == "__main__":
    unittest.main() 
