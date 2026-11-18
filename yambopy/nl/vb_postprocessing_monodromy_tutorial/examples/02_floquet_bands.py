# Floquet quasienergies of the propagated valence bands (and their m*hbar*omega
# replicas) compared with the Kohn-Sham bands along a high-symmetry path.
import matplotlib.pyplot as plt
from qepy.lattice import Path
from yambopy import FloquetAnalysis, FloquetBZ
#
fa = FloquetAnalysis(folder='.', calc='FLQ', save='SAVE')
bz = FloquetBZ('SAVE')            # full BZ of the (FixSymm) SAVE used in the RT run
#
path = Path([ [[  0.0,  0.0,  0.0],'$\\Gamma$'],
              [[  0.5,  0.0,  0.0],'M'],
              [[1./3.,1./3.,  0.0],'K'],
              [[  0.0,  0.0,  0.0],'$\\Gamma$']], [30,15,35])
#
data = fa.bands_on_path(path, npoints=40, bz=bz, ks_bands=range(1,9))   # raw numbers, if needed
fa.plot_bands(path, npoints=40, bz=bz, ks_bands=range(1,9), replicas=2, ylim=(-12,10),
              savefig='fig-Floquet_bands.png')
plt.show()
