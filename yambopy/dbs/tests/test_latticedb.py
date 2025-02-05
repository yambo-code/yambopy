# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
from qepy.lattice import Path
from yambopy.kpoints import get_path
from yambopy.dbs.latticedb import YamboLatticeDB
test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','gw_conv')

class TestYamboLatticeDB(unittest.TestCase):

    def test_yambolatticedb(self):
        """ test latticedb """

        #open latticedb
        filename = os.path.join(test_path,'SAVE/ns.db1')
        ydb = YamboLatticeDB.from_db_file(filename)

        #write json file
        ydb.write_json('lattice.json')        
        string1 = str(ydb)

        #read jsonfile
        ydb.from_json_file('lattice.json')
        string2 = str(ydb)
        
        assert string1 == string2

        p = Path([ [[1.0,1.0,1.0],'G'],
                   [[0.0,0.5,0.5],'X'],
                   [[0.0,0.0,0.0],'G'],
                   [[0.5,0.0,0.0],'L']], [20,20,20])
        bands_kpoints, bands_indexes, path_car = get_path(ydb.car_kpoints,ydb.rlat,None,p)

        print(ydb)

    def tearDown(self): 
        os.remove('lattice.json')

if __name__ == '__main__':
    unittest.main()
