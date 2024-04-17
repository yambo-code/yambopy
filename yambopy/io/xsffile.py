#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: RR
#
# This file is part of the yambopy project
#
import re
import numpy as np
from yambopy.lattice import replicate_red_kmesh, calculate_distances, get_path, car_red
from yambopy.dbs.latticedb import *

Bohr2Ang = 0.529177249

class YamboXsf():
    def __init__(self):
        self.cell_parameters = None
        self.cell_parameters_conv= None
        self.l_conventional = False
        self.dim = None
        self.natoms = None
        self.atom_positions = []
        self.grid_data = []
        self.data_array = np.array([])
        self.lattice = None
        self.grid_dim = np.array([])

    def set_cell_parameters(cls, cell_parameters):
        cls.cell_parameters = cell_parameters
    
    def set_cell_parameters_conv(cls, cell_parameters_conv):
        cls.cell_parameters_conv = cell_parameters_conv     

    def add_atom(cls, symbol, position):
        cls.atom_positions.append((symbol, position))

    def set_dim(cls,dim='3D'):
        cls.dim = dim
    
    def set_natoms(cls,natoms):
        cls.natoms = natoms        
        
    def set_lconventional(cls, l_conv):
        cls.l_conventional = l_conv
    
    def add_grid_data(cls, grid_name, block_dim, sub_grid_name, sub_grid_size, sub_grid_dim, sub_grid_origin,\
        sub_grid_vectors, data_array):
        cls.grid_data.append((grid_name, block_dim, sub_grid_name, sub_grid_size, sub_grid_dim, sub_grid_origin, \
            sub_grid_vectors, data_array))     

    def get_data_array(cls, data_array):
        cls.data_array = np.transpose(data_array,axes=(2,1,0)) # I want it in python order, to be checked if it's ok  

    def get_grid_dim(cls, grid_dim):
        cls.grid_dim = grid_dim

    def set_lattice(cls,lattice):
        cls.lattice = lattice

    def write_xsf(self, filename):
        with open(filename, 'w') as f:
            if (self.dim==3) : f.write("CRYSTAL\n")
            if (self.dim==2) : f.write("SLAB\n")
            if (self.dim==1) : f.write("POLYMER\n")
            if (self.dim==0) : f.write("MOLECULE\n")
            #PRIMVEC
            if np.any(self.cell_parameters):
                f.write("PRIMVEC\n")
                for vector in self.cell_parameters:
                    f.write("  ".join(str(x) for x in vector) + "\n")
            #CONVVEC
            if np.any(self.cell_parameters_conv):
                f.write("CONVVEC\n")
                for vector in self.cell_parameters_conv:
                    f.write("  ".join(str(x) for x in vector) + "\n")
            #ATOM_POSITIONS
            if (len(self.atom_positions) > 0):
                f.write("PRIMCOORD\n")
                f.write(f"{len(self.atom_positions)} 1\n")
                for symbol, position in self.atom_positions:
                    f.write(f"{symbol} {position[0]} {position[1]} {position[2]}\n")

        # wride grid_data for each sub_block
            if (len(self.grid_data)!=0):
                new_block_counter = 0 
                for grid_name, block_dim, sub_grid_name, sub_grid_size, sub_grid_dim, sub_grid_origin, sub_grid_vectors, sub_data_array in self.grid_data:
                    if(new_block_counter == len(sub_grid_name) or new_block_counter==0): f.write(f"BEGIN_BLOCK_DATAGRID_{block_dim}D\n")
                    f.write(f"\t{grid_name}\n")
                    f.write(f"\tBEGIN_DATAGRID_{sub_grid_dim}D\n")
                    f.write(f"".join("\t"+str(sub_grid_size[i]) for i in range(len(sub_grid_size)))+"\n")
                    f.write(f"".join("\t"+str(sub_grid_origin[i]) for i in range(len(sub_grid_origin)))+"\n")
                    for i in range(0,sub_grid_vectors.shape[0]):
                        f.write(f"".join("\t"+str(sub_grid_vectors[i][j]) for j in range(sub_grid_vectors.shape[1]))+"\n")
                    for i in range(sub_data_array.shape[0]):
                        for j in range(sub_data_array.shape[1]):
                            for k in range(sub_data_array.shape[2]):
                                #print(sub_data_array[i,j,k]) 
                                f.write(f"\t{sub_data_array[i,j,k]}\n")
                    f.write(f"\tEND_DATAGRID_{sub_grid_dim}D\n")
                    new_block_counter+=1
                    if(new_block_counter == len(sub_grid_name)or new_block_counter-1==0): f.write(f"END_BLOCK_DATAGRID_{block_dim}D")

    @staticmethod
    def read_xsf(filename):
        xsf_file = YamboXsf()

        with open(filename, 'r') as f:
            lines = [line.lstrip() for line in f.readlines()]

        #extract the dimension
        try: 
            if(lines.index("CRYSTAL\n")) : xsf_file.set_dim(3) 
            elif(lines.index("SLAB\n")) : xsf_file.set_dim(2)
            elif(lines.index("POLYMER\n")) : xsf_file.set_dim(1)
            elif(lines.index("MOLECULE\n")) : xsf_file.set_dim(0)
            else:
                ValueError("xsf file does not contain either CRYSTAL, SLAB, POLYMER, or MOLECULE")
        except ValueError as e:
            print(e)
        # Extract cell parameters
        cell_parameters = []
        cell_parameters_start = lines.index("PRIMVEC\n") + 1
        cell_parameters_end = cell_parameters_start + 3
        for line in lines[cell_parameters_start:cell_parameters_end]:
            cell_parameters.append([float(x) for x in line.split()])
        #Needed only if the convention cell is specified
        try:
            if(lines.index("CONVVEC\n")):
                xsf_file.set_lconventional(True)
                cell_parameters_conv=[]
                cell_parameters_conv_start = lines.index("CONVVEC\n")+1
                cell_parameters_conv_end = cell_parameters_conv_start+3
                for line in lines[cell_parameters_conv_start:cell_parameters_conv_end]:
                    cell_parameters_conv.append([float(x) for x in line.split()])    
                xsf_file.set_cell_parameters_conv(cell_parameters_conv)
            else:
                ValueError("conventional cell not present or not found.")
        except ValueError as e:       
            print(e) 
        # primitive cell parameters 
        xsf_file.set_cell_parameters(cell_parameters)
        # Extract atom positions
        atom_positions_start = lines.index("PRIMCOORD\n") + 1
        natoms = int(lines[atom_positions_start].split()[0]) # read the number of atoms
        atoms_positions_end = atom_positions_start +natoms +1
        for line in lines[atom_positions_start+1:atoms_positions_end]:
            symbol, x, y, z = line.split()
            position = [float(x), float(y), float(z)]
            xsf_file.add_atom(symbol, position)

        # Extract grid data
        block_start_idx = []
        block_dim = []
        
        for i,line in enumerate(lines): 
            if line.startswith(f"BEGIN_BLOCK_DATAGRID_"):
                block_start_idx.append(i)  
                block_dim.append(int(re.findall('[0-9]',line)[0])) 
                print(block_start_idx, block_dim)

        for istart_idx, start_idx in enumerate(block_start_idx):
            print(istart_idx,start_idx)
            block_end_idx = lines.index(f"END_BLOCK_DATAGRID_{block_dim[istart_idx]}D\n", start_idx) + 1
            block_lines = lines[start_idx:block_end_idx]
            #print(block_lines)

            subblock_start_idx = [i for i,line in enumerate(lines) if line.startswith(f"BEGIN_DATAGRID_{block_dim[istart_idx]}D\n")]
            print(subblock_start_idx)
            grid_name = block_lines[1].strip()
            for sub_start_idx in subblock_start_idx:
                sub_data_vectors=[]
                sub_block_end_idx = lines.index(f"END_DATAGRID_{block_dim[istart_idx]}D\n", sub_start_idx)+1
                sub_block_lines = lines[sub_start_idx:sub_block_end_idx]
                sub_block_name = re.split('_',sub_block_lines[0]) # extract name of the subblock
                sub_data_size = [int(x) for x in sub_block_lines[1].split()]
                sub_data_num = len(sub_data_size) # number of subdata for each block
                sub_data_origin = sub_block_lines[2].split()
                for j in range(3, 3+ sub_data_num):
                    sub_data_vectors.append([float(x) for x in sub_block_lines[j].split()])
                sub_data_values = []
                for j in range(3+sub_data_num,len(sub_block_lines)-1): # -1 becuase I do not want to include the last line
                    sub_data_values.append([float(x) for x in sub_block_lines[j].split()])
                sub_data_array = np.array(sub_data_values).reshape(np.flip(sub_data_size)) #reverse ordering cause fortran
                if(block_dim[istart_idx]==2): 
                    na = np.newaxis
                    sub_data_array = sub_data_array[:,:,na]
                xsf_file.get_data_array(data_array=sub_data_array)
                xsf_file.get_grid_dim(sub_data_size)
                xsf_file.add_grid_data(grid_name,block_dim=block_dim[istart_idx],sub_grid_name=sub_block_name, sub_grid_size=sub_data_size, \
                    sub_grid_dim=block_dim[istart_idx], sub_grid_origin=sub_data_origin, \
                    sub_grid_vectors=np.array(sub_data_vectors), data_array=sub_data_array)

        return xsf_file
    #instance of xsf file from lattice db class
    @staticmethod
    def from_latticedb(yambolattice):
        xsf_file = YamboXsf()
        xsf_file.set_dim(3) # I am not aware of Yambopy working with 2D 1D or 0D structures, for now
        xsf_file.set_lattice(yambolattice)
        xsf_file.set_cell_parameters(xsf_file.lattice.lat*Bohr2Ang)
        for iatom in range(0,len(xsf_file.lattice.car_atomic_positions)):
            xsf_file.add_atom(xsf_file.lattice.atomic_numbers[iatom],xsf_file.lattice.car_atomic_positions[iatom]*Bohr2Ang)
        return xsf_file
    
    def contribution_twolayers(self, zthreshold, c , fractional = True):
        data = self.data_array
        nGx, nGy, nGz = data.shape # G-vectors in real space
        # Now I need to know the range of G-vector in z that are < zthreshold or above
        # assume zthreshold is in fractionalcoordinates
        if (fractional):
            nGz_bottom = int(np.floor(nGz*(zthreshold)))
            bottom_layer = np.sum(data[:,:,:nGz_bottom])
            top_layer = np.sum(data[:,:,nGz_bottom:-1])
            return bottom_layer, top_layer
        else:
            zthreshold = zthreshold/c
            nGz_bottom = int(np.floor(nGz*(zthreshold)))
            print('ngz', nGz_bottom)
            bottom_layer = np.sum(data[:,:,:nGz_bottom])
            top_layer = np.sum(data[:,:,nGz_bottom:-1])
            return bottom_layer, top_layer
