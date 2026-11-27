# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
import numpy as np
from yambopy.bz.brillouin import BrillouinZone
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

        #k-points of the grid on a Brillouin-zone path
        bz = BrillouinZone(cell=ydb.lat, path_string='GXGL', npoints=60)
        bands_kpoints, bands_indexes, bands_distances = bz.find_collinear_grid(ydb)

        assert len(bands_kpoints) > 0
        assert np.all((0 <= bands_indexes) & (bands_indexes < ydb.ibz_nkpoints))
        assert np.all(np.diff(bands_distances) >= 0)
        assert np.allclose(bz.path_distances(bands_kpoints), bands_distances)

        print(ydb)

    def tearDown(self): 
        os.remove('lattice.json')

if __name__ == '__main__':
    unittest.main()
