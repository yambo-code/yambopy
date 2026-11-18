# Floquet analysis of one yambo_nl run (one driving frequency) with the
# direct monodromy method: quasienergies, Floquet-Kohn-Sham coefficients, checks.
import matplotlib.pyplot as plt
from yambopy import FloquetAnalysis
#
calc = 'FLQ'                      # folder with ndb.RT_V_bands, ndb.RT_V_bands_K_section, ndb.Nonlinear
fa = FloquetAnalysis(folder='.', calc=calc, save='SAVE')
print(fa)                         # summary and consistency checks
#
fa.report('quasienergies.dat')    # k, band, E_KS, xi, xi-E_KS, mode weight, edge weight, closure
fa.to_npz('floquet_FLQ.npz')      # everything needed to replot without the databases
#
kpt, band = 120, 4                # 1-based k-point index and propagated band
eta, d = fa.sambe(kpt, band)      # d~_kai(eta): array (2*N_F+1, n_basis)
xi = fa.quasienergies(gauge='ks') # (nk, nmodes) in eV, comparable with the KS energies
print(f'k={kpt} band={band}: E_KS = {fa.ks_energies()[kpt-1,1]:.6f} eV   xi = {xi[kpt-1,1]:.6f} eV')
#
fa.plot_sambe(kpt, band, savefig='fig-sambe.png')
fa.plot_periodic_part(kpt, band, basis_band=5, savefig='fig-periodic_part.png')
fa.plot_propagated_state(kpt, band, basis_band=5, savefig='fig-propagated_state.png')
fa.plot_checks(savefig='fig-checks.png')
plt.show()
