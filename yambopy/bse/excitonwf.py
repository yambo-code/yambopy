# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from __future__ import print_function, division
from yambopy import *
from itertools import product
from yambopy.plot import *
from yambopy.lattice import red_car

def jump_to(f,tag):
    """ Jump to a line in file
    """
    while True:
        line = f.readline()
        if tag in line:
            break
def v2str(v):
    return ("%12.8lf "*len(v))%tuple(v)

class YamboExcitonWaveFunctionXSF(object):
    """
    Class to read excitonic wavefunctions from yambo in the 3D xsf format
    """
    def __init__(self):
        self.cell = []
        self.atoms = []
        self.lattice = []
        self.initialized = False

    def read_file(self,filename):
        f = open(filename)
        jump_to(f,"PRIMVEC")
        self.lattice.append( list(map(float,f.readline().strip().split())) )
        self.lattice.append( list(map(float,f.readline().strip().split())) )
        self.lattice.append( list(map(float,f.readline().strip().split())) )

        jump_to(f,"PRIMCOORD")
        self.natoms = int(f.readline().split()[0])-1

        #read the hole position
        self.hole = list(map(float,f.readline().strip().split()))

        #read the atoms positions
        self.atoms = []
        for i in range(self.natoms):
            self.atoms.append( list(map(float,f.readline().strip().split())) )

        #get atypes
        self.atypes = np.unique([a[0] for a in self.atoms]).tolist()
        atypes_dict = dict([(a,n) for n,a in enumerate(self.atypes)])
        self.atoms = [ [atypes_dict[a[0]]]+a[1:] for a in self.atoms]

        jump_to(f,"BEGIN_DATAGRID_3D")
        self.nx, self.ny, self.nz = list(map(int, f.readline().strip().split()))
        f.readline() #ignore

        #read cell
        self.cell.append( list(map(float,f.readline().strip().split())) )
        self.cell.append( list(map(float,f.readline().strip().split())) )
        self.cell.append( list(map(float,f.readline().strip().split())) )

        #read data
        self.datagrid = np.zeros([self.nz,self.ny,self.nx])
        for k,j,i in product(list(range(self.nz)),list(range(self.ny)),list(range(self.nx))):
            self.datagrid[k,j,i] = float(f.readline())
        self.initialized = True

    def plot_slice_x(self,n):
        """ plot a slice of the 3d grid
        """
        plt.imshow(self.datagrid[:,:,n])
        plt.show()

    def plot_slice_z(self,n):
        """ plot a slice of the 3d grid
        """
        plt.imshow(self.datagrid[n,:,:])
        plt.show()

    def write_xsf(self,filename):
        f = open(filename,'w')
        #structure
        f.write('CRYSTAL\n')
        f.write('PRIMVEC\n')
        for vlat in self.lattice:
            f.write(v2str(vlat)+'\n')
        f.write('PRIMCOORD\n')
        f.write('%d 1\n'%len(self.atoms))
        f.write(v2str(self.hole)+'\n')
        for atom in self.atoms:
            f.write("%d "%self.atypes[atom[0]]+v2str(atom[1:])+'\n')
        #datagrid
        f.write('BEGIN_BLOCK_DATAGRID_3D\n')
        f.write('excitonwf\n')
        f.write('BEGIN_DATAGRID_3D\n')
        f.write('%d %d %d\n'%(self.nx,self.ny,self.nz))
        f.write('0.00 0.00 0.00\n')
        for vlat in self.lattice:
            f.write(v2str(vlat)+'\n')
        for x in self.datagrid.flatten():
            f.write('%lf\n'%x)
        f.write('END_DATAGRID_3D\n')
        f.write('END_BLOCK_DATAGRID_3D\n')
        f.close()

    def center(self):
        """ Center the atoms and excitonic wavefunction in the unit cell
        """
        self.norm = [0,0,0]
        for i in range(3):
            self.norm[i] = np.linalg.norm(self.lattice[i])

        print(self.norm)
        print(self.nx, self.ny, self.nz)

        # find the average position in each direction (geometric center)
        pos = np.zeros([3])
        for atom in self.atoms:
            pos += atom[1:]
        pos = pos/len(self.atoms)

        # center the atoms around that position using the center of the unit cell
        displecement = pos
        for atom in self.atoms:
            x,y,z = atompos + atom[1:]
            atom[1:] = np.array([x%self.nx, y%self.ny, z%self.nz])

        # center the wavefunction around that position
        new_data = np.zeros(data.shape)

    def get_data(self):
        return { "datagrid": self.datagrid.flatten().tolist(),
                 "lattice": self.lattice,
                 "atoms": self.atoms,
                 "atypes": self.atypes,
                 "hole": self.hole,
                 "nx": self.nx,
                 "ny": self.ny,
                 "nz": self.nz }

    def write_json(self):
        """ Write as a json file
        """
        JsonDumper(self.get_data(),"datagrid.json")

    def read_json_file(self,filename):
        f = open(filename,"r")
        data = json.load(f)
        f.close()

        self.read_json(data)

    def read_json(self,data):
        """ Write as a json file
        """
        self.datagrid = data["datagrid"]
        self.lattice = data["lattice"]
        self.atoms = data["atoms"]
        self.atypes = data["atypes"]
        self.nx = data["nx"]
        self.ny = data["ny"]
        self.nz = data["nz"]
        self.natoms = len(self.atoms)

        self.initialized = True

    def __str__(self):
        s = ""
        s += "lattice:\n"
        for i in range(3):
            s += ("%12.8lf "*3)%tuple(self.lattice[i])+"\n"
        s += "natoms:\n"
        s += "atoms:\n"
        for i in range(self.natoms):
            s += ("%3d "+"%12.8lf "*3)%tuple(self.atoms[i])+"\n"
        s += "atypes:\n"
        for n,a in enumerate(self.atypes):
            s += "%3d %3d\n"%(n,a)
        s += "nx: %d\n"%self.nx
        s += "ny: %d\n"%self.ny
        s += "nz: %d\n"%self.nz
        return s
