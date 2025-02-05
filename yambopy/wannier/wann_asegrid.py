import sisl
import spglib
from yambopy.lattice import red_car
from yambopy.units import bohr2ang
import numpy as np
from yambopy.wannier.wann_io import NNKP
from yambopy.wannier.wann_kpoints import KPointGenerator
from ase.dft.kpoints import monkhorst_pack
from yambopy.units import *

class ase_Monkhorst_Pack(KPointGenerator):
    def __init__(self, grid_shape,latdb):
        super().__init__()
        self.grid_shape = grid_shape
        self.latdb = latdb
        self.rlat = self.latdb.rlat*2*np.pi 
    
    def generate(self):
        # Use ASE's monkhorst_pack to generate the k-points
        self.k = monkhorst_pack(self.grid_shape)
        self.nkpoints = len(self.k)
        self.weights = 1/(self.nkpoints)
        self.red_kpoints = self.k
        self.car_kpoints = red_car(self.k, self.rlat)*ang2bohr # result in Bohr
        print(f"Generated {self.nkpoints} k-points using ASE.")