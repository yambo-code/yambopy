# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import numpy as np
import unittest
import os
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.electronsdb import YamboElectronsDB
from qepy.lattice import Path

test_path = os.path.join(os.path.dirname(__file__),'..','..','data','refs','bse')

class TestYamboExcitonDB(unittest.TestCase):

    def test_yamboexcitondb(self):

        #define path in reduced coordinates
        npoints = 20
        path = Path([ [[0.0, 0.0, 0.0],'$\\Gamma$'],
                      [[0.5, 0.0, 0.0],'M'],
                      [[1./3,1./3,0.0],'K'],
                      [[0.0, 0.0, 0.0],'$\\Gamma$']],
                      [int(npoints*2),int(npoints),int(np.sqrt(5)*npoints)])

        #load databases
        lat  = YamboLatticeDB.from_db_file(os.path.join(test_path,'SAVE','ns.db1'))
        exc  = YamboExcitonDB.from_db_file(lat,folder=os.path.join(test_path,'yambo'))
        electrons = YamboElectronsDB(lat,save=os.path.join(test_path,'SAVE')) 

        #show output
        sort_e, sort_i = exc.get_sorted()
        exc.get_degenerate(1)

        #write exc_I.dat and exc_E.dat
        exc.write_sorted('exc')

        #calculate chi
        w,chi = exc.get_chi()

        #load reference
        eps_filename = os.path.join(test_path,'o-yambo.eps_q1_diago_bse')
        w_ref,chi_imag_ref,chi_real_ref = np.loadtxt(eps_filename,unpack=True)[:3]

        #get amplitude and phases
        kpoints, amplitude, phase = exc.get_amplitudes_phases((1,0,))

        #exciton_bandstructure
        exc.plot_exciton_bs(electrons, path, (1,2,), args_plot={'c':'g'},space='bands',show=False)

    def tearDown(self):
        if os.path.isfile('exc_I.dat'): os.remove('exc_I.dat')
        if os.path.isfile('exc_E.dat'): os.remove('exc_E.dat')


if __name__ == '__main__':
    unittest.main()
