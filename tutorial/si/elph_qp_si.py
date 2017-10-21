#############################################################################
#
# Author: Alejandro Molina-Sanchez
# Run electron-phonon calculations using Yambo 
#
# Calculations are done inside the folder elphon 
#
##############################################################################
#from __future__ import print_function
from yambopy import *
from qepy    import *
from numpy import loadtxt,ones
import argparse

#parse options
parser = argparse.ArgumentParser(description='Convergence test of the colomb cutoff')
parser.add_argument('-r' ,'--run',      action="store_true", help='Run the calculation')
parser.add_argument('-p' ,'--plot',     action="store_true", help='Run the analysis')
parser.add_argument('-a' ,'--advanced', action="store_true", help='Run the analysis')
args = parser.parse_args()

yambo_ph  = 'yambo_ph'
ypp_ph    = 'ypp_ph'

folder_ya = 'elphon'

temperature = [0,500]

# A. Run QPs El-Ph correction
if args.run:

  # Check the existence of el-ph matrix elements
  if not os.path.isdir('%s/SAVE'%folder_ya):
    print('Electron-phonon matrix elemenst are missing...')
    print('Run script elph_pw_si.py')
    exit()

  # Generatio of yambo input-file 
  yqp = YamboIn('yambo_ph -g n -c p',folder=folder_ya,filename='yambo.in')
  ysf = YamboIn('yambo_ph -g g -c p',folder=folder_ya,filename='yambo.in')
  yqp.arguments.append('ExtendOut')

  # Calculation changing temperature
  for T in temperature:
    # QP correction
    yqp['BoseTemp'] = [ T, 'K' ] 
    # Spectral Functions
    ysf['BoseTemp'] = [ T, 'K' ] 
    ysf['QPkrange'] = [ 1, 1, 4, 5 ] 
    ysf['GEnRnge']  = [ [-5, 5], 'eV' ] 
    ysf['GEnSteps'] = 500 
    yqp.write('%s/qp-%dK.in'% (folder_ya,T))
    ysf.write('%s/sf-%dK.in'% (folder_ya,T))
    os.system('cd %s; yambo_ph -F qp-%dK.in -J qp-%dK'%(folder_ya,T,T))
    os.system('cd %s; yambo_ph -F sf-%dK.in -J sf-%dK'%(folder_ya,T,T))

# B. Analysis

  #pack the files in .json files
if args.plot:
  pack_files_in_folder(folder_ya)
  #plot the results using yambo analyser
  ya = YamboAnalyser()
  print('plot QPs corrections')
  ya.plot_qp_correction('qp')
  #path = [[[0.5,   0,   0],'L'],
  #        [[  0,   0,   0],'$\Gamma$'],
  #        [[  0, 0.5, 0.5],'X'],
  #        [[1.0, 1.0, 1.0],'$\Gamma$']]
  #ya.plot_gw_path('qp',path,cols=(3,))
  print('plot Spectral Functions')
  ya.plot_spectral_function('band')

# C. Advance options

if args.advanced:
  # Run
  yqp.arguments.append('WRgFsq')
