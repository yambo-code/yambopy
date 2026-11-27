from qepy import PwXML
from yambopy import BrillouinZone

# k-points map
npoints = 50
intervals = [npoints,npoints,npoints,npoints,npoints]
extra_points = {'G': [0.0, 0.0, 0.0], 'H': [0.0, 0.0, 1.0], 'N': [1/2, 0.0, 1/2], 'P': [1/2, 1/2, 1/2]}
path = 'GHNGPN'

bz = BrillouinZone(
    path_string=path,
    extra_points=extra_points,
    intervals=intervals
)

# Class PwXML. QE database reading
xml = PwXML(prefix='pw', path='bands')
xml.plot_eigen(bz.get_indices())
