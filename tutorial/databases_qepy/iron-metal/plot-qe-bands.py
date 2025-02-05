from qepy import *

# k-points map
npoints = 50
path_kpoints = Path([ [[0.0, 0.0, 0.0 ],'G'],
                      [[0.0, 0.0, 1.0 ],'H'],
                      [[1./2,0.0,1./2.],'N'],
                      [[0.0, 0.0, 0.0 ],'G'],
                      [[1./2, 1./2, 1./2 ],'P'],
                      [[1./2,0.0,1./2. ],'N']], [npoints,npoints,npoints,npoints,npoints])

# Class PwXML. QE database reading
xml = PwXML(prefix='pw',path='bands/t0')

# Class PwXML. QE database reading
xml.plot_eigen(path_kpoints)
