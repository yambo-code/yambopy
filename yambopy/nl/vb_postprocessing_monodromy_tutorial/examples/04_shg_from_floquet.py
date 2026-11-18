# Second-harmonic generation reconstructed from the Floquet states
# (Berry-phase polarization of the reconstructed states) and comparison with the
# real-time result. Requires ndb.Overlap in every FLQ folder.
import glob
import numpy as np
import matplotlib.pyplot as plt
from yambopy import FloquetAnalysis, FloquetBZ, YamboNLDB, floquet_harmonic_spectrum
from yambopy.nl.sin_analysis import Xn_from_sine
#
# --- a single frequency, step by step
bz = FloquetBZ('SAVE')
fa = FloquetAnalysis(folder='.', calc='FLQ', save='SAVE')
nldb = YamboNLDB(calc='FLQ')
fa.setup_berry(bz=bz, field_direction=nldb.Efield[0]['versor'])
print(fa.validate_berry(nldb))                  # Berry polarization of the stored states vs ndb.Nonlinear
chi2, info = fa.harmonic_response(nldb, order=2)
print(f"{info['energy_eV']:.4f} eV  chi2 (x,y,z) = {chi2}  [cm/statV]")
#
# --- the whole Floquet frequency grid (one folder per frequency)
calcs = sorted(glob.glob('FLQ_[0-9][0-9][0-9]'))
E_fl, X_fl = floquet_harmonic_spectrum(calcs, save='SAVE', order=2)
#
# --- real-time reference (multi-frequency run in folder NL)
rt = YamboNLDB(calc='NL')
ana = Xn_from_sine(rt, X_order=4)
freqs, X_rt = ana.output_analysis(ana.perform_analysis(), to_file=False)
#
plt.plot(freqs*27.2113834, abs(X_rt[2,:,1]), label='real time')
plt.plot(E_fl, abs(X_fl[:,1]), 'o', label='Floquet')
plt.xlabel('Driving energy (eV)'); plt.ylabel(r'$|\chi^{(2)}_y|$ (cm/statV)'); plt.legend()
plt.savefig('fig-SHG_RT_vs_Floquet.png')
plt.show()
