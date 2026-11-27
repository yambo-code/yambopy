from qepy import PwXML
from yambopy import BrillouinZone
from math import sqrt

# k-points map
npoints = 50
intervals = [int(npoints*2), int(npoints), int(sqrt(5)*npoints)]
extra_points = {'G': [0.0, 0.0, 0.0], 'M': [0.5, 0.0, 0.0], 'K': [1/3, 1/3, 0.0]}
path = 'GMKG'

bz = BrillouinZone(
    path_string=path,
    extra_points=extra_points,
    intervals=intervals
)

# Class PwXML. QE database reading
xml = PwXML(prefix='bn', path='bands')
xml.plot_eigen(bz.get_indices())

# Alternative option with more plot control
## Matplotlib options 
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

## Class PwXML. QE database reading
xml.plot_eigen_ax(ax,bz.get_indices())
plt.show()
