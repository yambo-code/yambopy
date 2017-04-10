#
#
#
#
import os
import re
import math
import numpy as np
from qepy import PwIn
from itertools import product
import copy
from math import *

class supercell():
    """A class to generate custom supercells from a quantum espresso input file
    """

    def __init__(self,qe_input,R,type='diagonal'):
        """ 
        qe_input: a PwIn() instance of an input file in the unit cell (uc)
        R: 
            if type='diagonal' R is a list of the repetitions of the uc in the cartesian directions (integers)
            else R contains the fractional coordinates of the q-point to be folded at Gamma in a nondiagonal supercell
        """
        self.qe_input = qe_input
        self.R        = R
        self.latvec   = qe_input.cell_parameters
        self.basis    = int(qe_input.system['nat'])
        self.atoms    = qe_input.atoms
        self.sup_size = R[0]*R[1]*R[2]

        if type=='diagonal': new_latvec,new_atoms = self.build_diagonal()
        else:   self.build_nondiagonal()
        
        #PwIn() object that can be printed, written to file, etc.
        self.qe = self.write(new_atoms,new_latvec)

    def build_diagonal(self):
        latvec     = self.latvec
        R          = self.R
        new_latvec = np.array([R[0]*latvec[0],R[1]*latvec[1],R[2]*latvec[2]])
        atoms      = np.array([atom[1] for atom in self.atoms])
        #new_atoms[cell][basis][direction]
        new_atoms      = np.array([atoms for n in range(self.sup_size)])
        for nz,ny,nx,b in product(range(R[2]),range(R[1]),range(R[0]),range(self.basis)): 
            cell=nx+ny*R[0]+nz*R[0]*R[1]
            new_atoms[cell,b,0]=(new_atoms[cell,b,0]+nx)/R[0]
            new_atoms[cell,b,1]=(new_atoms[cell,b,1]+ny)/R[1]
            new_atoms[cell,b,2]=(new_atoms[cell,b,2]+nz)/R[2]
        #new_atoms[super_basis][directions]$
        new_atoms=new_atoms.reshape(self.basis*self.sup_size,3)
        return new_latvec, new_atoms
     
    def build_nondiagonal(self):
        """Nondiagonal supercell, based on [REFERENCE]
        """
        exit()    
    
    def atoms_input(self, new_atoms):
        """ Put the atomic element labels in the right order
        """
        positions_input = new_atoms.tolist()
        elements_input  = [[self.qe_input.atoms[i][0] for i in range(self.basis)] for j in range(self.sup_size)]
        elements_input  = [ item for sublist in elements_input for item in sublist ]
        atoms_input     = [[elements_input[i], positions_input[i]] for i in range(self.sup_size*self.basis)]
        return atoms_input
    
    def write(self,new_atoms,new_latvec):
        R = self.R
        qe = self.qe_input
        #Just a suggestion for the new kpoint mesh 
        new_kpoints = [ceil(qe.kpoints[0]/R[0]), ceil(qe.kpoints[1]/R[1]), ceil(qe.kpoints[2]/R[2])]
        qe_s = copy.deepcopy(qe)
        qe_s.atoms = self.atoms_input(new_atoms)
        qe_s.control['prefix'] = qe.control['prefix'][:-1]+"_s'"
        if 'celldm(1)' in qe.system: qe_s.system['celldm(1)'] = R[0]*float(qe.system['celldm(1)'])
        if 'celldm(2)' in qe.system: qe_s.system['celldm(2)'] = R[1]*float(qe.system['celldm(2)'])/float(qe_s.system['celldm(1)'])
        if 'celldm(3)' in qe.system: qe_s.system['celldm(3)'] = R[2]*float(qe.system['celldm(3)'])/float(qe_s.system['celldm(1)'])
        if qe.cell_parameters: qe_s.cell_parameters = new_latvec
        #Just a suggestion for the new bands
        if 'nbnd' in qe.system: qe_s.system['nbnd'] = self.sup_size*qe.system['nbnd']
        qe_s.system['nat'] = self.basis*self.sup_size
        qe_s.kpoints = new_kpoints
        return qe_s
