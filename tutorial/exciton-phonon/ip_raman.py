"""
IP one-phonon Stokes Raman with yambopy.

The driver does all the band-window bookkeeping:
  - loads lattice, electrons, dipoles, ndb.elph
  - BZ-expands the dipoles via YamboDipolesDB(expand=True) and reorders
    them into the LetzElPhC k-grid
  - slices electron energies, dipoles and gkkp to the user-chosen
    `raman_bands` window (Fortran 1-indexed, inclusive)
  - passes ready-to-use arrays to ip_resonant_raman_tensor_oneph

The Raman function itself takes only the energies, dipoles, gkkp,
cell volume and broadening and evaluates the formula.

Assumes yambo wrote velocity-gauge dipoles (DIP_v in ndb.dipoles).
"""

import numpy as np

from yambopy import (YamboLatticeDB, YamboElectronsDB,
                     YamboDipolesDB, LetzElphElectronPhononDB)
from yambopy.exciton_phonon.excph_resonant_raman import ip_resonant_raman_tensor_oneph


# ---- User inputs ----------------------------------------------------------
folder      = '.'
SAVE_dir    = folder + '/SAVE'
elph_path   = folder + '/ndb.elph'
dipole_path = folder + '/dipoles/ndb.dipoles'

omega_range = (1.0, 6.0, 501)         # laser energies (eV): min, max, npts
broading    = 0.10                    # Lorentzian broadening (eV)
ph_fre_th   = 5.0                     # acoustic-mode cutoff (cm^-1)
raman_bands = (44, 56)                # 1-indexed Fortran inclusive


# ---- Databases ------------------------------------------------------------
lattice   = YamboLatticeDB.from_db_file(filename='%s/ns.db1' % SAVE_dir)
electrons = YamboElectronsDB.from_db_file(folder=SAVE_dir,
                                          filename='ns.db1', Expand=True)
dipdb     = YamboDipolesDB.from_db_file(lattice,
                                        filename=dipole_path,
                                        dip_type='v',
                                        project=False,
                                        expand=True)
elph      = LetzElphElectronPhononDB(elph_path, read_all=True)


# ---- Yambo BZ -> elph BZ permutation --------------------------------------
diff = elph.kpoints[:, None, :] - lattice.red_kpoints[None, :, :]
diff = diff - np.round(diff)
perm = np.argmin(np.linalg.norm(diff, axis=-1), axis=1)


# ---- Resolve band range (1-indexed inclusive, Fortran) -------------------
b_in_F, b_out_F = raman_bands
nb              = b_out_F - b_in_F + 1
n_val           = int(electrons.nbandsv) - (b_in_F - 1)
nc              = nb - n_val
cell_vol        = float(lattice.lat_vol)


# ---- eph_g(q=0) sliced to the raman window --------------------------------
elph_b_in_F = int(min(elph.bands))
off         = b_in_F - elph_b_in_F

iq0   = 0
g_q0  = elph.gkkp[iq0, :, :, 0, :, :]                       # (nk, nmodes, nb_i, nb_f)
g_q0  = np.swapaxes(g_q0, -1, -2)
eph_g = np.transpose(g_q0, (1, 0, 2, 3)) / 2.0              # (nmodes, nk, nb_eph, nb_eph), Ha
eph_g = eph_g[:, :, off : off + nb, off : off + nb]         # (nmodes, nk, nb, nb)

ph_eV = elph.ph_energies[iq0]                               # (nmodes,), eV


# ---- KS energies in elph k-order, sliced to the raman window --------------
el_eV = electrons.eigenvalues[0][perm, b_in_F - 1 : b_out_F]   # (nk, nb)


# ---- Velocity dipoles: cv block, reorder, slice to raman cv block --------
nv_dip       = int(dipdb.nbandsv)
nc_dip       = int(dipdb.nbandsc)
elec_dipoles = dipdb.dipoles[:, :, nv_dip : nv_dip + nc_dip, :nv_dip]   # (nk_y, 3, nc_dip, nv_dip)
elec_dipoles = elec_dipoles[perm]                                       # -> elph k-order
elec_dipoles = elec_dipoles[:, :, :nc, -n_val:]                         # raman cv block
elec_dipoles = np.transpose(elec_dipoles, (1, 0, 2, 3))                 # (3, nk, nc, n_val)


# ---- Raman calculation ----------------------------------------------------
laser_eV = np.linspace(omega_range[0], omega_range[1], int(omega_range[2]))

R = ip_resonant_raman_tensor_oneph(
    laser_energies    = laser_eV,
    ph_energies       = ph_eV,
    el_energies       = el_eV,
    elec_dipoles      = elec_dipoles,
    eph_g             = eph_g,
    cell_vol          = cell_vol,
    broad             = broading,
    ph_freq_threshold = ph_fre_th,
)


# ---- Save -----------------------------------------------------------------
raman_inte = np.sum(np.abs(R)**2, axis=(1, 2, 3))

np.savetxt('raman_intensity.dat',
           np.column_stack([laser_eV, raman_inte]),
           header='laser_energy(eV)   raman_intensity(arb.u.)')

np.savez('Raman_tensors.npz',
         laser_energies_eV = laser_eV,
         ph_freq_cm        = ph_eV * 8065.54,
         raman_tensor      = R)

print('Wrote raman_intensity.dat and Raman_tensors.npz')
