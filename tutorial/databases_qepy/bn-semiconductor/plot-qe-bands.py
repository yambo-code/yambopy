from qepy import *
from math import sqrt

# k-points map
npoints = 50
path_kpoints = Path([ [[0.0, 0.0, 0.0],'$\Gamma$'],
                      [[0.5, 0.0, 0.0],'M'],
                      [[1./3,1./3,0.0],'K'],
                      [[0.0, 0.0, 0.0],'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)])

atom_1 = [range(16)]
atom_2 = [range(16,32)]

# Class PwXML. QE database reading
xml = PwXML(prefix='bn',path='bands')

# Class PwXML. QE database reading
xml.plot_eigen(path_kpoints)
