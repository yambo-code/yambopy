from __future__ import print_function
#
#
#
#
import os
import re
import math
import numpy as np
from qepy import PwIn
from yambopy.lattice import *
from itertools import product
import copy
from math import *
import fractions as frc

class supercell():
    """A class to generate custom supercells from a quantum espresso input file
    """

    def __init__(self,qe_input,R,mode='diagonal',units='fractional',write=True):
        """ 
        qe_input: a PwIn() instance of an input file in the unit cell (uc)
        R: 
        if default: R is a list of the repetitions of the uc in the cartesian directions (integers)
        else: R contains the fractional coordinates of the q-point to be folded at Gamma in a nondiagonal supercell like [[m1,m2,m3],[n1,n2,n3]]
        units: atomic positions in fractional(/angstrom/bohr)
        """
        self.qe_input = qe_input
        self.latvec   = np.array(qe_input.cell_parameters)
        self.basis    = int(qe_input.system['nat'])
        self.atoms    = qe_input.atoms
        self.b2a      = 0.529177
        #Case of nondiagonal supercell
        if mode!='diagonal':
            self.Q = np.array(R)
            print('Nondiagonal supercell')
            if (qe_input.kpoints % self.Q[1] != 0).any():
                print('ERROR: You must set a unit cell k-point mesh where%s\
       Nx,Ny,Nz are multiples of %d,%d,%d, respectively.'%('\n',self.Q[1,0],self.Q[1,1],self.Q[1,2])) 
                exit()
            self.R, self.new_latvec = self.find_nondiagonal()
        #Case of diagonal supercell    
        else: 
            self.R = R
            self.sup_size = R[0]*R[1]*R[2]
            self.new_latvec = np.array([self.latvec[i]*R[i] for i in range(3)])

        self.sup_size = self.R[0]*self.R[1]*self.R[2]
        new_atoms = self.build_supercell()        
        if write:
            #PwIn() object that can be printed, written to file, etc.
            self.qe = self.write(new_atoms,mode)

    def lattice_constants(self,vec):
        return [np.linalg.norm(vec[0]),np.linalg.norm(vec[1]),np.linalg.norm(vec[2])]

    def build_supercell(self):
        latvec     = self.latvec
        R          = self.R
        atoms      = np.array([atom[1] for atom in self.atoms])
        atoms      = red_car(atoms,latvec) 
        #new_atoms[cell][basis][direction]
        new_atoms      = np.array([atoms for n in range(self.sup_size)])
        for nz,ny,nx,b in product(range(R[2]),range(R[1]),range(R[0]),range(self.basis)): 
            cell=nx+ny*R[0]+nz*R[0]*R[1]
            new_atoms[cell,b]=new_atoms[cell,b] +nx*latvec[0] +ny*latvec[1] +nz*latvec[2]
        #new_atoms[super_basis][directions]$
        new_atoms=new_atoms.reshape(self.basis*self.sup_size,3)
        new_atoms=car_red(new_atoms,self.new_latvec)
        return new_atoms

    def find_integers(self,nums,g23,g12,g31,g123):
        """Compute integers for off-diagonal supercell matrix elements 
           Called by find_nondiagonal()
        """
        #Compute p (it's a modulo equation)
        if g23 == 1: p = 0
        else:
            for i in range(1,g23):
                if (nums[1]+i*nums[2]) % g23 == 0:
                    p=i
                    break
         #Compute q
        g12_r = g12/g123
        g23_r = g23/g123
        g31_r = g31/g123
        if g12_r == 1: q = 0
        else:
            for i in range(1,g12_r):
                if (g23_r*nums[0]+i*g31_r*nums[1]) % g12_r == 0:
                    q=i
                    break
        #Compute r
        gg_r = g31*g23/g123
        z = g23*nums[0]/g12+g31*q*nums[1]/g12
        if gg_r == 1: r = 0
        else:
            for i in range(1,gg_r):
                if (z+i*nums[2]) % gg_r == 0:
                    r=i
                    break
        return p,q,r 
 
    def find_nondiagonal(self):
        """Nondiagonal supercell, based on [Phys. Rev. B 92, 184301]
        """
        Q = self.Q
        #Take care of components already at Gamma
        Q[1,np.where(Q[0]==0)]=1
        #Shift the q-point into the positive quadrant of the reciprocal unit cell
        Q[0,np.where(Q[0]<0)]+=Q[1,np.where(Q[0]<0)]
        #GCDs of Q[1] (in the logical order of the derivation)
        g23  = frc.gcd(Q[1,1],Q[1,2])
        g12  = frc.gcd(Q[1,0],Q[1,1])
        g31  = frc.gcd(Q[1,2],Q[1,0])
        g123 = frc.gcd(Q[1,0],frc.gcd(Q[1,1],Q[1,2]))
        #Integers needed to solve the supercell matrix equation    
        p,q,r = self.find_integers(Q[0],g23,g12,g31,g123)            
        #Matrix elements (in order of derivation) and supercell matrix
        S_33 =        Q[1,2]
        S_22 =        Q[1,1]/g23
        S_23 =      p*Q[1,2]/g23
        S_11 =   g123*Q[1,0]/(g12*g31)
        S_12 = q*g123*Q[1,1]/(g12*g23)
        S_13 = r*g123*Q[1,2]/(g31*g23)
        self.S = np.array([[S_11,S_12,S_13],[0,S_22,S_23],[0,0,S_33]])
        #New lattice vectors and actual supercell size
        new_latvec = np.einsum('ij,jx->ix',self.S,self.latvec)
        R          = [self.S[0,0],self.S[1,1],self.S[2,2]]
        print(self.S)
        return R, new_latvec

    def reciprocal(self,mode):
        """Function to compute reciprocal lattice
        """
        #Unit cell
        repvec = rec_lat(self.latvec)
        alat=np.array(self.lattice_constants(self.latvec))
        self.repvec = 2.*np.pi*np.multiply(1./alat,repvec)
        #Supercell
        if mode=='diagonal': self.new_repvec = np.array([self.repvec[i]/float(R[i]) for i in range(3)])
        else: 
            self.S_inv_T = np.linalg.inv(self.S).T
            self.new_repvec = np.einsum('ij,jx->ix',self.S_inv_T,self.repvec)

    def atoms_input(self, new_atoms):
        """ Put the atomic element labels in the right order
        """
        positions_input = new_atoms.tolist()
        elements_input  = [[self.qe_input.atoms[i][0] for i in range(self.basis)] for j in range(self.sup_size)]
        elements_input  = [ item for sublist in elements_input for item in sublist ]
        atoms_input     = [[elements_input[i], positions_input[i]] for i in range(self.sup_size*self.basis)]
        return atoms_input

    def posint(self,value):
        return abs(int(round(value)))

    def write(self,new_atoms,mode):
        R = self.R
        new_latvec = self.new_latvec
        alat = self.lattice_constants(new_latvec)
        qe = self.qe_input
        if mode=='diagonal':
            #A suggestion for a consistent new kpoint mesh 
            new_kpoints = [ceil(qe.kpoints[0]/R[0]), ceil(qe.kpoints[1]/R[1]), ceil(qe.kpoints[2]/R[2])]
        else:
            #The compulsory new kpoint mesh - (sub)multiples of it are also fine but not consistent
            self.reciprocal('nondiagonal')
            new_kpoints = np.dot(self.S_inv_T,np.array(qe.kpoints))
            new_kpoints = [self.posint(new_kpoints[0]),self.posint(new_kpoints[1]),self.posint(new_kpoints[2])]
        qe_s = copy.deepcopy(qe)
        qe_s.atoms = self.atoms_input(new_atoms)
        qe_s.control['prefix'] = qe.control['prefix'][:-1]+"_s'"
        if 'celldm(1)' in qe.system:
            qe_s.system['celldm(1)'] = alat[0]
            qe_s.system['celldm(2)'] = alat[1]/alat[0]
            qe_s.system['celldm(3)'] = alat[2]/alat[0]
        if not (mode=='diagonal'and R[0]==R[1]==R[2]): qe_s.system['ibrav']=0
        qe_s.cell_units = 'bohr'
        qe_s.cell_parameters = new_latvec
        #Just a suggestion for the new bands
        if 'nbnd' in qe.system: qe_s.system['nbnd'] = self.sup_size*int(qe.system['nbnd'])
        qe_s.system['nat'] = self.basis*self.sup_size
        qe_s.kpoints = new_kpoints
        return qe_s
