# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
from yambopy.dbs.qpdb import YamboQPDB
test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','gw')

class TestYamboQPDB(unittest.TestCase):

    def test_yamboqpdb(self):

        #open qpdb
        qpdb = YamboQPDB.from_db(folder=test_path)

        #get qp dbb
        eigenvalues_qp, eingenvalues_dft, lifetimes, z = qpdb.get_qps()

        #plot bs
        qpdb.plot_bs(show=False)

        print(qpdb)

if __name__ == '__main__':
    unittest.main()
