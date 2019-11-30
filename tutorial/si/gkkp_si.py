#
# Author: Alejandro Molina-Sanchez 
# Obtain the electron-phonon matrix elements for Real-time simulations
# 
#
from __future__ import print_function
from yambopy import *
from qepy import *

ph = 'ph.x'
folder = 'elphon'

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

os.system('mkdir -p %s'%folder)

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
phin['nq1'], phin['nq2'], phin['nq3'] = 2, 2, 2

phin.write('%s/02ph.in'%folder)      # Potential calculation
phin['trans']           = '.false.'
phin['electron_phonon'] = "'yambo'"
phin.write('%s/04elph.in'%folder)    # Electron-phonon calculation

# A. Generation s.dbph_# Files

# 1. Self-consistent data
os.system('cp -r scf/si.save %s'%folder)
# 2. Potential dVscf 
os.system('cd %s; %s < 02ph.in   | tee 02.out'%(folder,ph))
# 3. Non-self consistent data
os.system('cp -r nscf/si.save %s/.'%folder)
# 4. Electron-phonon matrix elements
os.system('cd %s; %s < 04elph.in | tee 04.out'%(folder,ph))

# B. Generation of the gkkp fragments 

# 1. Database in Yambo
os.system('cd nscf/si.save; p2y -O ../../%s/ELPH'%folder)
os.system('cp %s/elph_dir/* %s/ELPH'%(folder,folder))
# 2. Setup yambo
y  = YamboIn.from_runlevel('yambo_rt -i -V all -Q',folder='%s/ELPH'%folder)
y.arguments.append('BSEscatt')
y.write('%s/ELPH/yambo.in'%folder)
os.system('cd %s/ELPH ; yambo_rt -F yambo.in'%folder)
# 3. Expansion gkkp matrix elements
yp = YamboIn.from_runlevel('ypp_ph -g',folder='%s/ELPH'%folder,filename='ypp.in')
yp.arguments.append('GkkpExpand')
yp.write('%s/ELPH/ypp.in'%folder)
os.system('cd %s/ELPH ; ypp_ph -F ypp.in'%folder)
# 4. Moving files to rt
os.system('mkdir -p rt; cd rt ; mkdir -p GKKP')
os.system('mv %s/ELPH/SAVE/ndb.elph* rt/GKKP'%folder)
print('Files ready in folder rt')
