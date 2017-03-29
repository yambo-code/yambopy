from __future__ import print_function
import sys
from yambopy     import *
from qepy        import *
from schedulerpy import *
import argparse


# GW plot
#pack_files_in_folder('gw')
ygw = YamboAnalyser('gw')
path = [[[0,   0,   0],'$\Gamma$'],
        [[0.5, 0,   0],'M'],
        [[0.3333,0.3333, 0.0],'K'],
        [[0.0, 0.0, 0.0],'$\Gamma$']]
#ygw.plot_gw_path('qp',path, cols=(lambda x: x[2]+x[3],))

# BSE plot
#ybse = YamboOut('bse')
#ybse.pack()
#ybse.plot(tag='eps', cols=(2,))

exc = YamboExcitonWeight(path='bse',filename='bse/o-yambo.exc_weights_at_try',save='bse/SAVE')

exc.plot_exciton_band(path)
