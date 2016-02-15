#
# Author: Alejandro Molina-Sanchez 
# Obtain the electron-phonon matrix elements for Real-time simulations
# 
#
from __future__ import print_function
from yambopy import *
from qepy import *

ph = 'ph.x'

#check if the scf cycle is present
if os.path.isdir('scf/si.save'):
    print('scf calculation found!')
else:
    print('scf calculation not found!')
    exit()

#check if the nscf cycle is present
if os.path.isdir('nscf/si.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit()

# Create a work directory

os.system('mkdir -p work_elph')

# Input files for the electron-phonon calculation

phin = PhIn()
phin['prefix']          = "'si'"
phin['tr2_ph']          = 1.0e-8
phin['fildyn']          = "'si.dyn'"
phin['fildvscf']        = "'dvscf'"
phin['iverbosity']      = 1
phin['ldisp']           = '.true.'
phin['trans']           = '.true.'
phin['electron_phonon'] =  "'dvscf'" 
phin['nq1']             = 2
phin['nq2']             = 2
phin['nq3']             = 2

phin.write('work_elph/02ph.in')      # Potential calculation
phin['trans']           = '.false.'
phin['electron_phonon'] = "'yambo'"
phin.write('work_elph/04elph.in')    # Electron-phonon calculation

# A. Generation s.dbph_# Files

# 1. Self-consistent data
os.system('cp -r scf/si.save work_elph')
# 2. Potential dVscf 
os.system('cd work_elph; %s < 02ph.in   | tee 02.out'%ph)
# 3. Non-self consistent data
os.system('cp -r nscf/si.save work_elph/.')
# 4. Electron-phonon matrix elements
os.system('cd work_elph; %s < 04elph.in | tee 04.out'%ph)

# B. Generation of the gkkp fragments 

# 1. Database in Yambo
os.system('cd nscf/si.save; p2y -O ../../work_elph/ELPH')
os.system('cp work_elph/elph_dir/* work_elph/ELPH ; cd work_elph/ELPH')
# 2. Setup yambo
y  = YamboIn('yambo_rt -i -V all -Q',folder='work_elph/ELPH')
y.arguments.append('BSEscatt')
y.write('work_elph/ELPH/yambo.in')
os.system('cd work_elph/ELPH ; yambo_rt -F yambo.in')
# 3. Expansion gkkp matrix elements
yp = YamboIn('ypp_ph -g',folder='work_elph/ELPH',filename='ypp.in')
yp.arguments.append('GkkpExpand')
yp.write('work_elph/ELPH/ypp.in')
os.system('cd work_elph/ELPH ; ypp_ph -F ypp.in')
# 4. Moving files to FixSymm
os.system('mkdir -p FixSymm; cd FixSymm ; mkdir -p GKKP')
os.system('mv work_elph/ELPH/SAVE/ndb.elph* FixSymm/GKKP')
print('Files ready in folder FixSymm')
