from netCDF4 import Dataset
from yambopy import YamboVbandsDB, VbPP
calc='fl_sample'
mydb=YamboVbandsDB(calc=calc)
print(mydb)
#
kpt,band = 6,4
vbs = VbPP(mydb,kpt,band)
print(vbs)
#
optimized_evecs = vbs.find_qe()
print(f'QE= {optimized_evecs.FL_qe} eV with accuracy {optimized_evecs.nr_acc} eV after {optimized_evecs.nr_it} iterations - error in periodicity = {optimized_evecs.err}')
optimized_evecs.plot_realtime(band_to_plot=1,t_step=0.004)
optimized_evecs.plot_floquet()
optimized_evecs.output()
