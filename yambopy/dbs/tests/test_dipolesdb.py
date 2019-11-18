# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
from yambopy.dbs.dipolesdb import YamboDipolesDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.electronsdb import YamboElectronsDB

test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','ip','SAVE')

class TestYamboDipolesDB(unittest.TestCase):

    def test_yambodipolesdb(self):

        # read lattice
        lat = YamboLatticeDB.from_db_file(os.path.join(test_path,'ns.db1'))

        #read electrons
        electrons = YamboElectronsDB(lat,save=test_path) 
        print(electrons)

        #read dipoles
        # P  -> the velocity matrix elements
        # iR -> from the velocity matrix elements using r = [H,r]
        dipoles = YamboDipolesDB(lat,save=test_path,dip_type='P')

        #calculate epsilon
        dipoles.ip_eps2(electrons)

        print(dipoles)

if __name__ == '__main__':
    unittest.main()
