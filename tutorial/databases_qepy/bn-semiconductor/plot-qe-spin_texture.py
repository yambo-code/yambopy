from qepy import *
from yambopy import *

# Figure features
fig = plt.figure(figsize=(7,7))
ax = fig.add_axes( [ 0.15, 0.1, 0.8, 0.8 ])

# DB definitions
prefix = 'bn'
path = '.'
folder_spin = 'spin_bands'
prefix_spin = 'spin'

# Initializing Spin_texture class
Spin_tex = Spin_texture(prefix=prefix,path=path,folder_spin=folder_spin,prefix_spin=prefix_spin)

# Load spin projection values (i = 1,2,3 = x,y,z)
Spin_x = Spin_tex.load_spin_data(1)

print(Spin_x)

# Plot spin_texture (mode = "raw", "interpolated", "arrow")
Spin_tex.plot_spin_texture(ax,limfactor=0.8,mode="raw",nband=5)

plt.show()



