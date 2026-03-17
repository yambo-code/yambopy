from netCDF4 import Dataset
from yambopy import YamboVbandsDB, YamboLatticeDB, findallk_qe, plot_2D_kdist
from qepy.lattice import Path
import numpy as np
import matplotlib.pyplot as plt
#
calc='fl_sample'
mydb=YamboVbandsDB(calc=calc)
fl_eignvectors, fl_quasienergy, ks_eigenvalues = findallk_qe(mydb,report_file='results_eval.dat')
#
lat = YamboLatticeDB.from_db_file(filename='SAVE/ns.db1',Expand=False)
floq_mode, bnd_1, bnd_2 = 2,4,5
floq_vbnd = bnd_1 - 1 
floq_cbnd = bnd_2 -  mydb.basis_index[0]
#
data = np.zeros([mydb.n_kpts])
data = np.absolute(fl_eignvectors[:,floq_vbnd,floq_mode,floq_cbnd])
#
kdist = plot_2D_kdist(data,lat,nspin=-1,plt_cbar=True,shift_BZ=False)
