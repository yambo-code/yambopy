from yambopy import *

# 1. pack files of convergence GW calculations
# important: all calculations must finish
pack_files_in_folder('bse_calc')

# 2. Read json files
ya = YamboAnalyser('bse_calc')

# 3. Plot eps-BSE
ya.plot_bse('eps_q1',cols=(2,))
