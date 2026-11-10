import numpy as np
import matplotlib.pyplot as plt
from qepy import PwXML,Unfolding
from yambopy import BrillouinZone

# Matplotlib options
fig = plt.figure(figsize=(7,5))
ax = fig.add_axes( [ 0.1, 0.1, 0.4, 0.8 ])

# k-points map
npoints = 30

# Path of the supercell bands calculation, M-G-K-M with npoints per segment.
# Only the tick positions and labels are needed, so no lattice is given.
bz_sc = BrillouinZone(path_string='MGKM', intervals=[npoints, npoints, npoints])

# Prefix and path to read data
prefix_pc = 'BN_pc'
prefix_sc = 'BN_sc_2x2x2'

path_pc = 'Unfolding/BN_pc'
path_sc = 'Unfolding/BN_sc_2x2x2'

# Reading QE database
pc = PwXML(prefix=prefix_pc,path=path_pc)
sc = PwXML(prefix=prefix_sc,path=path_sc)

# Computing projections of sc values into pc
fold = Unfolding(prefix_pc=prefix_pc,path_pc=path_pc,prefix_sc=prefix_sc,path_sc=path_sc,spin="noncol", band_min = 0, sc_rotated=False,compute_projections=True)

# Plotting unfolded sc band structure
fold.plot_eigen_ax(ax,path=bz_sc.get_indices(),ylim = (-4.0,10.0))

plt.show()

