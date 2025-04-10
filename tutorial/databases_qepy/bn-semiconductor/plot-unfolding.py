import numpy as np
from qepy *
import matplotlib.pyplot as plt

# Matplotlib options
fig = plt.figure(figsize=(7,5))
ax = fig.add_axes( [ 0.1, 0.1, 0.4, 0.8 ])

# k-points map
npoints = 30

path_kpoints_pc = Path([ [[0.0,  0.5, 0.0],'M'],
                         [[0.0,  0.0, 0.0],'G'],
                         [[1.0/3.0,  1.0/3.0, 0.0],'K'],
                         [[0.0,  0.5, 0.0],'G']], [int(npoints),int(npoints),int(npoints)])

path_kpoints_sc = Path([ [[0.0,  1.0, 0.0],'M'],
                         [[0.0,  0.0, 0.0],'G'],
                         [[2.0/3.0,  2.0/3.0, 0.0],'K'],
                         [[0.0,  1.0, 0.0],'G']], [int(npoints),int(npoints),int(npoints)])

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
fold.plot_eigen_ax(ax,path=path_kpoints_sc,ylim = (-4.0,10.0))

plt.show()

