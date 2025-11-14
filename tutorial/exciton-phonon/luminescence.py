"""
This script computes the phonon-assisted luminescence spectra for indirect materials
"""
import numpy as np
import matplotlib.pyplot as plt
from yambopy.exciton_phonon.excph_luminescence import exc_ph_luminescence
from yambopy.exciton_phonon.excph_input_data import exc_ph_get_inputs

path = '3D_hBN'
# Path to BSE calculation (Lout--> response is Lfull)
bsepath =  f'{path}/bse_Lfull' # ndb.BS_diago_Q* databases are needed
# Path to BSE calculation for optically active exciton (Lin --> response is Lbar)
bseBARpath =  f'{path}/bse_Lbar'  # ndb.BS_diago_Q1 database is needed
# Path to electron-phonon calculation
elphpath = path # ndb.elph is needed
# Path to unprojected dipoles matrix elements (optional)
dipolespath = bsepath # ndb.dipoles is needed (optional)
# Path to lattice and k-space info
savepath = f'{path}/SAVE' # ns.db1 database is needed

bands_range=[6,10] # 2 valence, 2 conduction bands
phonons_range=[0,12] # All phonons
nexc_out = 12 # 12 excitonic states at each momentum (Lout)
nexc_in  = 12 # 12 excitonic states at Q=0 (Lin)
T_ph  = 10 # Lattice temperature
T_exc = 10 # Effective excitonic temperature

emin=4.4      # Energy range and plot details (in eV)
emax=4.7
estep=0.0002
broad = 0.005 # Broadening parameter for peak width (in eV)

# We calculate and load all the inputs:
# * Exciton-phonon matrix elements
# * Excitonic dipole matrix elements
# * Exciton energies
# * Phonon energies
# * We specify bse_path2=bseBARpath meaning we use Lbar calculation for Q=0 excitons
input_data = exc_ph_get_inputs(savepath,elphpath,bsepath,\
                               bse_path2=bseBARpath,dipoles_path=dipolespath,\
                               nexc_in=12,nexc_out=12,\
                               bands_range=[6,10],phonons_range=[0,12])

ph_energies, exc_energies, exc_energies_in, G, exc_dipoles = input_data

# We calculate the luminescence spectrum including the input data from before 
w,PL = exc_ph_luminescence(T_ph,ph_energies,exc_energies,exc_dipoles,G,\
                           exc_energies_in=exc_energies_in,exc_temp=T_exc,\
                           nexc_out=nexc_out,nexc_in=nexc_in,emin=emin,emax=emax,\
                           estep=estep,broad=broad) #,ph_channels='e')

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_xlim(emin,emax)
ax.set_ylim(0,np.max(PL)*1.1)
ax.get_yaxis().set_visible(False)

ax.plot(w, PL, '-',c='red', label="AA' hBN luminescence")

plt.legend()
plt.savefig('hBN_luminescence.png')
plt.show()

