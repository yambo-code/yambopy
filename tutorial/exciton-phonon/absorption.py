"""
This script computes the phonon-assisted optical absorption (dielectric function)
for indirect materials, including exciton-phonon satellites (Eqs. 1-3).
"""
import numpy as np
import matplotlib.pyplot as plt
from yambopy import YamboLatticeDB
from yambopy.exciton_phonon.excph_absorption import exc_ph_absorption
from yambopy.exciton_phonon.excph_input_data import exc_ph_get_inputs

path = '3D_hBN'
# Path to BSE calculation (Lout--> response is Lfull)
bsepath =  f'{path}/bse_Lfull' # ndb.BS_diago_Q* databases are needed
# Path to BSE calculation for optically active exciton (Lin --> response is Lbar)
bseBARpath =  f'{path}/bse_Lbar'  # ndb.BS_diago_Q1 database is needed
# Path to electron-phonon calculation
elphpath = path # ndb.elph is needed
# Path to unprojected dipoles matrix elements
dipolespath = bsepath # ndb.dipoles is needed
# Path to lattice and k-space info
savepath = f'{path}/SAVE' # ns.db1 database is needed

phonons_range=[0,12] # All phonons
nexc_out = 12 # 12 excitonic states at each momentum (Lout)
nexc_in  = 12 # 12 excitonic states at Q=0 (Lin)
T_ph  = 10    # Lattice temperature (no excitonic temperature: absorption starts from the ground state)

estep   = 0.0002
broad   = 0.005  # eta: HWHM of the Lorentzians in eV
broad_0 = 0.001  # delta: imaginary part in the satellite denominators in eV

renorm_direct     = True  # (1-R) factor on the direct BSE peaks
renorm_satellites = True  # (1-R) weight of the intermediate exciton in the satellites
pol_dir = None            # None: average over x,y,z; e.g. [1,0,0] for in-plane light

# Unit cell volume (bohr^3) for the absolute prefactor 4*pi*s/Omega
lattice = YamboLatticeDB.from_db_file(filename=f'{savepath}/ns.db1')
cell_volume = lattice.lat_vol

# We calculate and load all the inputs (same as for luminescence):
# * Exciton-phonon matrix elements
# * Excitonic dipole matrix elements
# * Exciton energies
# * Phonon energies
# * We specify bse_path2=bseBARpath: Lbar at Q=0 gives the macroscopic dielectric function
input_data = exc_ph_get_inputs(savepath,elphpath,bsepath,\
                               bse_path2=bseBARpath,dipoles_path=dipolespath,\
                               nexc_in=nexc_in,nexc_out=nexc_out,\
                               phonons_range=phonons_range)

ph_energies, exc_energies, exc_energies_in, G, exc_dipoles = input_data

# Energy range: from below the lowest indirect satellite to above the direct excitons (eV)
emin = np.min(exc_energies) - np.max(ph_energies) - 0.1
emax = np.max(exc_energies_in) + np.max(ph_energies) + 0.2

# Phonon-assisted absorption (mode-resolved satellites)
w, eps, deps_modes, R = exc_ph_absorption(T_ph,ph_energies,exc_energies,exc_dipoles,G,cell_volume,\
                                          exc_energies_in=exc_energies_in,\
                                          nexc_out=nexc_out,nexc_in=nexc_in,emin=emin,emax=emax,\
                                          estep=estep,broad=broad,broad_0=broad_0,pol_dir=pol_dir,\
                                          renorm_direct=renorm_direct,renorm_satellites=renorm_satellites,\
                                          mode_resolved=True,return_renorm=True) #,ph_channels='e')
deps     = np.sum(deps_modes,axis=1)  # total phonon-assisted correction
eps_dir  = eps - deps                 # direct BSE term (renormalised if renorm_direct)

# Bare BSE (no exciton-phonon coupling) for comparison
_, eps_bse, _ = exc_ph_absorption(T_ph,ph_energies,exc_energies,exc_dipoles,np.zeros_like(G),cell_volume,\
                                  exc_energies_in=exc_energies_in,nexc_out=nexc_out,nexc_in=nexc_in,\
                                  emin=emin,emax=emax,estep=estep,broad=broad,pol_dir=pol_dir)

# Save dielectric functions
print("Printing phonon assisted absorption in 'absorption.dat' file")
np.savetxt('absorption.dat',np.c_[w,eps.imag,eps.real,deps.imag,eps_dir.imag,eps_bse.imag],\
           header='E[eV]  Im(eps_tot)  Re(eps_tot)  Im(Deps_ph)  Im(eps_direct)  Im(eps_BSE_bare)')
np.savetxt('absorption_modes.dat',np.c_[w,deps_modes.imag],\
           header='E[eV]  Im(Deps_ph) for each phonon mode')
np.savetxt('renormalization.dat',R.T,header='R per optical exciton (columns: phonon channels)')

fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(6,7))

# Linear scale: full spectrum
ax1.plot(w, eps_bse.imag, '--', c='gray',  label='BSE')
ax1.plot(w, eps.imag,     '-',  c='blue',  label='BSE + exc-ph')
ax1.fill_between(w, deps.imag, color='red', alpha=0.3, label='phonon satellites')
ax1.set_ylabel(r'Im[$\varepsilon(\omega)$]')
ax1.set_ylim(0,np.max(eps.imag)*1.1)
ax1.legend()

# Log scale: low-energy tail (phonon-assisted onset) and mode decomposition
ax2.semilogy(w, eps_bse.imag, '--', c='gray', label='BSE')
ax2.semilogy(w, eps.imag,     '-',  c='blue', label='BSE + exc-ph')
for im in range(deps_modes.shape[1]):
    if np.max(deps_modes[:,im].imag) > 1e-3*np.max(deps.imag):
        ax2.semilogy(w, deps_modes[:,im].imag, lw=0.8, label=f'mode {phonons_range[0]+im+1}')
ax2.set_ylim(np.max(deps.imag)*1e-3, np.max(eps.imag)*2)
ax2.set_xlim(emin,emax)
ax2.set_xlabel('Energy (eV)')
ax2.set_ylabel(r'Im[$\varepsilon(\omega)$]')
ax2.legend(fontsize=7,ncol=2)

plt.tight_layout()
plt.savefig('hBN_absorption.png')
plt.show()
