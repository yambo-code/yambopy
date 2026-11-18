# Weight of the Floquet sidebands on the conduction bands in the Brillouin zone:
#   W_eta(k) = sum_a sum_{i in conduction} |d~_kai(eta)|^2
# eta = +-1: one-photon (linear) component, eta = +-2: two-photon (SHG) component.
import matplotlib.pyplot as plt
from qepy.lattice import Path
from yambopy import FloquetAnalysis, FloquetBZ
#
fa = FloquetAnalysis(folder='.', calc='FLQ', save='SAVE')
bz = FloquetBZ('SAVE')
path = Path([ [[0.0,0.0,0.0],'$\\Gamma$'], [[0.5,0.0,0.0],'M'],
              [[1./3.,1./3.,0.0],'K'], [[0.0,0.0,0.0],'$\\Gamma$']], [30,15,35])
#
conduction = [5, 6]
for eta in (1, 2):
    fa.plot_sideband_bz(eta, basis_bands=conduction, bz=bz, path=path,
                        savefig=f'fig-sideband_eta{eta}_BZ.png')
#
w2 = fa.sideband_weight(2, basis_bands=conduction)        # (nk,) values on the IBZ
w2_full = bz.expand(w2)                                   # same on the full BZ
print('k-points with the largest eta=2 weight:', (w2.argsort()[::-1][:5]+1).tolist())
plt.show()
