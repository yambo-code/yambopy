
# This script checks convergence for all k-points

from yambopy import *
import numpy as np
import matplotlib.pyplot as plt

# pack files of convergence GW calculations
pack_files_in_folder('gw_conv')

# Start Analyser
ya = YamboAnalyser('gw_conv')

# Plot of all the k-points converging one parameter
ya.plot_gw_all_kpoints_convergence(tag='EXX')

ya.plot_gw_all_kpoints_convergence(tag='Bnds')

ya.plot_gw_all_kpoints_convergence(tag='NGsBlk')

ya.plot_gw_all_kpoints_convergence(tag='GbndRnge')
