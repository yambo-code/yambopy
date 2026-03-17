from netCDF4 import Dataset
from yambopy import YamboVbandsDB, YamboLatticeDB, findallk_qe, calc_rho, plot_2D_kdist,get_qebands_path, get_qebands_interpolate
from qepy.lattice import Path
import numpy as np
import matplotlib.pyplot as plt
#
calc='fl_sample'  # name of the directory containing the relevant DBs
mydb=YamboVbandsDB(calc=calc)
fl_eignvectors, fl_quasienergy, ks_eigenvalues = findallk_qe(mydb,report_file='results_eval.dat')
#
npoints = 10  # number of divisions
path = Path([ [[  0.0,  0.0,  0.0],'G'],  # path along the BZ
              [[  0.5,  0.0,  0.0],'M'],
              [[1./3.,1./3.,  0.0],'K'],
              [[  0.0,  0.0,  0.0],'G']], [int(npoints*2),int(npoints),int(np.sqrt(5)*npoints)] )
#
lat = YamboLatticeDB.from_db_file(filename='SAVE/ns.db1',Expand=False)
ks_bs, fl_bs = get_qebands_interpolate(lat,path,fl_quasienergy, ks_eigenvalues)
fig = plt.figure(figsize=(4,5))
ax = fig.add_axes( [ 0.20, 0.20, 0.70, 0.70 ])
ks_bs.plot_ax(ax,legend=True,c_bands='r',ylim=(-6.5,0.0),label='KS')
fl_bs.plot_ax(ax,legend=True,c_bands='b',linestyle= 'dashed',ylim=(-6.5,0.0),label='FL')
plt.savefig('QE_band_interp.pdf')
plt.close()
