# Copyright (C) 2017 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
import unittest
import os
from yambopy.dbs.latticedb import YamboLatticeDB

test_path = os.path.join('..','..','..','tests','reference','gw_conv')

class TestYamboLatticeDB(unittest.TestCase):

    def test_yambolatticedb(self):
        """ test latticedb """

        #open latticedb
        filename = os.path.join(test_path,'SAVE/ns.db1')
        ydb = YamboLatticeDB.from_db_file(filename)

        #write json file
        ydb.write_json('lattice.json')        

        #read jsonfile
        ydb.from_json_file('lattice.json')

    def tearDown(self): 
        os.remove('lattice.json')

if __name__ == '__main__':
    unittest.main()
