from __future__ import print_function
#############################################################################
#
# Author: Alejandro Molina-Sanchez
# Run electron-phonon calculations using Yambo 
#
# This script run QE to obtain the electron-phonon matrix elements to
# be used by Yambo
#
# Calculations are done inside the folder elphon 
#
##############################################################################
#from __future__ import print_function
from yambopy import *
from qepy    import *
from numpy import loadtxt,ones

yambo_ph  = 'yambo_ph'
ypp_ph    = 'ypp_ph'
ph        = 'ph.x'

folder_ph = 'work'
folder_ya = 'elphon'
nq        = 10 

# check if the database is present
if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
if os.path.isdir('nscf/si.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit() 

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    os.system('cd nscf/si.save; p2y')
    os.system('cd nscf/si.save; yambo')
    os.system('mv nscf/si.save/SAVE database')

os.system('mkdir -p %s'%folder_ph)

#if not os.path.  ramdom
def qpoint_generation(nqpoints):
  if os.path.isfile('%s/qlist.dat'%folder_ph):
      print('List of q-points already exists')
      return loadtxt('%s/qlist.dat'%folder_ph)
  else:
      print('creating random list of q-points')
      ypp = YamboIn('ypp_ph -k r',folder='database',filename='ypp.in')
      ypp['cooOut'] = "alat"
      ypp['BZ_random_Nk']= nqpoints
      ypp.arguments.append('NoWeights')
      ypp.write('database/ypp.in')
      os.system('cd database; %s -F ypp.in'%ypp_ph)
      q        = ones((nqpoints,4))
      q[:,:-1] = loadtxt('database/o.random_k_pts')
      q[0]     = [ 0.0, 0.0, 0.0, 1.0]
      os.system('rm database/ypp.in') 
      os.system('rm database/o.random_k_pts')
      f = open('%s/qlist.dat'%folder_ph,'w')
      for iq in q:
        f.write('%12.8f '*4%tuple(iq) + '\n')
      f.close()
      print('List of q-random points in alat units')
      return q
  
# Number of random qpoints

qlist = qpoint_generation(nq)

phin = PhIn()
phin['prefix']          = "'si'"
phin['tr2_ph']          = 1.0e-5
phin['fildyn']          = "'si.dyn'"
phin['fildvscf']        = "'dvscf'"
phin['iverbosity']      = 1
phin['ldisp']           = '.false.'
phin['trans']           = '.true.'
phin['qplot']           = '.true.'
phin['electron_phonon'] =  "'dvscf'" 
phin.qpoints            = qlist 

# A. Writing files for ph.x 

phin.write('%s/02ph.in'%folder_ph)      # Potential calculation
phin['trans']           = '.false.'
phin['electron_phonon'] = "'yambo'"
phin.write('%s/04elph.in'%folder_ph)    # Electron-phonon calculation

# B. Generation s.dbph_# Files (check if they exist)

qcount = 0
for iq in range(1,nq+1):
  print('s.dbph_%06d'%iq)
  if os.path.isfile('%s/elph_dir/s.dbph_%06d'%(folder_ph,iq)):
    qcount += 1
if qcount == nq:
  print('El-ph matrix already calculated')
else:
  # 1. Self-consistent data
  os.system('cp -r scf/si.save %s'%folder_ph)
  # 2. Potential dVscf 
  os.system('cd %s; %s < 02ph.in   | tee 02.out'%(folder_ph,ph))
  # 3. Non-self consistent data
  os.system('cp -r nscf/si.save %s/.'%folder_ph)
  # 4. Electron-phonon matrix elements
  os.system('cd %s; %s < 04elph.in | tee 04.out'%(folder_ph,ph))

# C. Generation El-Ph DBs for yambo

os.system('mkdir -p %s'%folder_ya)
if not os.path.isdir('%s/SAVE'%folder_ya):
  os.system('mkdir %s/SAVE'%folder_ya)
  os.system('cp -r database/SAVE/ns.* %s/SAVE/.'%folder_ya)
  os.system('cp %s/elph_dir/s.dbph* %s/.'       %(folder_ph,folder_ya))

if not os.path.isfile('%s/SAVE/ndb.elph_gkkp'%folder_ya):
  yin = YamboIn('yambo_ph -i -V kpt',folder=folder_ya,filename='yambo.in')
  yin.arguments.append('MinusQ')
  yin.write('%s/init.in'%folder_ya)
  os.system('cd %s ; yambo_ph -F init.in'%folder_ya)

  ypp = YamboIn('ypp_ph -g',folder=folder_ya,filename='ypp.in')
  os.system('cd %s; ypp_ph'%folder_ya)
