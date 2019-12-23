from yambopy import *

# 1. pack files of convergence GW calculations
# important: all calculations must finish
pack_files_in_folder('bse_conv')

# 2. Read json files
ya = YamboAnalyser('bse_conv')

# Plot BSE spectra for each parameter
ya.plot_bse(('eps_q1','FFTGvecs'),cols=(2,),png_file=True)

ya.plot_bse(('eps_q1','NGsBlkXs'),cols=(2,),png_file=True)

ya.plot_bse(('eps_q1','BndsRnXs'),cols=(2,),png_file=True)
