"""
Tutorial for YamboSpinDB and band magnetic quantities.

Read spin, orbital moment and Berry curvature from a stock ndb.dipoles.

EDIT the paths below to point to the yambo SAVE folder and the dipoles
database. The dipoles database must be generated with DipBandsAll.
"""

import numpy as np

from yambopy import YamboDipolesDB, YamboElectronsDB, YamboLatticeDB
from yambopy.dbs.berry_curvature import compute_electronic_berry_curvature
from yambopy.dbs.orbital_magnetic_moment import (
    compute_electronic_magnetic_moment,
    compute_electronic_orbital_moment,
)
from yambopy.dbs.spindb import YamboSpinDB


# Simple defaults for running this script from the calculation folder.
save_path = "SAVE"
dipoles_file = "stock_dipoles_db/ndb.dipoles"
bands_range = [23, 30]
kpoint = 37       # one-based IBZ k-point index
band = 26         # one-based band index


# Load the original Yambopy databases.
lattice = YamboLatticeDB.from_db_file(filename=save_path + "/ns.db1")
electrons = YamboElectronsDB.from_db_file(folder=save_path)
dipoles = YamboDipolesDB.from_db_file(
    lattice, filename=dipoles_file, dip_type="v",
    expand=False, project=False)

# Read spin from the same stock dipole database.
spin = YamboSpinDB.from_db_file(dipoles_file, bands_range=bands_range)

print("Loaded databases")
print("  dipoles.dipoles_full shape :", np.shape(dipoles.dipoles_full))
print("  spin.spin_full shape       :", np.shape(spin.spin_full))
print("  electron energies shape    :", np.shape(electrons.eigenvalues_ibz))
print("  dipole band window         :", dipoles.bands_range)
print("  spin band window           :", spin.bands_range)

# Compute band-level quantities.
orbital = compute_electronic_orbital_moment(
    dipoles, electrons, bands_range=bands_range)
magnetic = compute_electronic_magnetic_moment(
    dipoles, electrons, filename=dipoles_file,
    bands_range=bands_range, spindb=spin)
berry = compute_electronic_berry_curvature(
    dipoles, electrons, bands_range=bands_range)

ik = kpoint - 1
ib = band - bands_range[0]

sigma_z = spin.spin_full[ik, 2, ib, ib].real
sz_over_hbar = 0.5*sigma_z
lz_over_hbar = orbital["L_over_hbar"][2, ik, 0, ib]
mu_orb = orbital["mu_orb_over_muB"][2, ik, 0, ib]
mu_spin = magnetic["spin_moment_over_muB"][2, ik, 0, ib, ib].real
mu_total = magnetic["total_moment_over_muB"][2, ik, 0, ib, ib].real
omega_z = berry["Omega_ang2"][2, ik, 0, ib]

print("")
print("Selected quantity")
print("  k-point                  : %d (one-based IBZ)" % kpoint)
print("  band                     : %d (one-based)" % band)
print("  sigma_z                  : % .8f" % sigma_z)
print("  S_z/hbar = sigma_z/2     : % .8f" % sz_over_hbar)
print("  L_z/hbar                 : % .8f" % lz_over_hbar)
print("  mu_orb,z/mu_B = -L_z/hbar: % .8f" % mu_orb)
print("  mu_spin,z/mu_B           : % .8f" % mu_spin)
print("  mu_total,z/mu_B          : % .8f" % mu_total)
print("  Omega_z                  : % .8f Ang^2" % omega_z)

if not orbital["valid"][ik, 0, ib]:
    print("")
    print("Warning: this band/k-point is near-degenerate within degen_tol.")
    print("The diagonal orbital and Berry values are diagnostics only.")
