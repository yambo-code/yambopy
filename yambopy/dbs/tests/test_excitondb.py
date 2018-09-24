# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.electronsdb import YamboElectronsDB

test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','bse')

class TestYamboExcitonDB(unittest.TestCase):

    def test_yamboexcitondb(self):

        #define path in reduced coordinates
        path = [ [0.0, 0.0, 0.0],
                 [0.5, 0.0, 0.0],
                 [1./3,1./3,0.0],
                 [0.0, 0.0, 0.0]]

        #load databases
        lat  = YamboLatticeDB.from_db_file(os.path.join(test_path,'SAVE','ns.db1'))
        exc  = YamboExcitonDB(lat,path=os.path.join(test_path,'yambo'))
        electrons = YamboElectronsDB(lat,save=os.path.join(test_path,'SAVE')) 

        #show output
        print(exc)
        exc.get_sorted()
        exc.get_degenerate(1)

        #write exc_I.dat and exc_E.dat
        exc.write_sorted('exc')

        #get amplitude and phases
        kpoints, amplitude, phase = exc.get_amplitudes_phases((1,0,))

        #exciton_bandstructure
        #import matplotlib.pyplot as plt
        #ax = plt.gca()
        #exc.plot_exciton_bs(ax, lat, electrons, path, (1,2,), args_plot={'c':'g'},space='bands')

    def tearDown(self):
        os.remove('exc_I.dat')
        os.remove('exc_E.dat')


if __name__ == '__main__':
    unittest.main()
