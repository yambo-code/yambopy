#
# 1st version by Fulvio Paleari
# This file is part of yambopy
#
import os
import re
import numpy as np
from qepy import PwIn
from yambopy.lattice import *
from itertools import product
import copy
from math import *
# Dimensional constants for reference (actually only b2a is really used in the code for now)
cm1_2_Tera=0.0299793         # Conversion from cm-1 to THz with 2pi factor included
Tera=1.e12                   
b2a =0.529177                
hbar=6.5821e-16              # Planck's constant (eV*s)
kb=8.6173e-5                 # Boltzmann's constant (eV/K)
Mp=1.0073                    # Proton mass (reference, u)
cMp=Mp*1.660539*6.241509e-29 # Conversion of Mp in eV*\AA^{-2}*s^2
#
## ISSUES TO FIX ##
"""(iii) Small issue in nondiagonal supercell matrices for certain q-vectors
   (iv)  In nondiagonal supercells: the suggested kpoint mesh must be CONSISTENT and UNIFORM --> to fix
"""
## These two functions read phonons from qe output
def read_frequencies(modes_file,units='Tera'):
    """Read phonon frequencies from QE output phonon modes file
    """
    Omega=[]
    with open (modes_file) as fp:
        for line in fp:
            if line.strip()[0:4]=='freq':
                w=re.findall(r"[-+]?\d*\.\d+|d+", line)
                Omega.append(w)
    Omega = np.float_(Omega)
    if units=='Tera': Omega= Omega[:,0]
    else:             Omega= Omega[:,1]
    return Omega

def read_eig(modes_file,basis):
    """ Read phonon modes from QE output file
    """
    modes=3*basis
    eig = np.genfromtxt(modes_file, autostrip=True, comments='freq', skip_header=4, skip_footer=1, usecols=(1,2,3,4,5,6))
    # Reshape data from quantum espresso output
    eig = np.array([eig[:,0]+1j*eig[:,1], eig[:,2]+1j*eig[:,3], eig[:,4]+1j*eig[:,5]])
    eig = eig.T
    eig = np.reshape(eig, (modes,basis,3))
    # eig[mode][atom][direction]
    return eig

## Start of the proper supercell class
class Supercell():
    """
        A class to generate custom supercells from a quantum espresso input file.
        The supercell may be "non-diagonal" according to [Lloyd-Williams and Monserrat, Phys. Rev. B 92, 184301, 2015].
        The atoms in the supercell may also be displaced along phonon modes if provided.
        
        Input arguments:
       
            - qe_input: a PwIn() instance of an input file in the unit cell (uc)
            - Optional: QE-DFPT phonon modes output
        
        How it works:
            
            - Call self.d_sup(R) to build diagonal supercell
                -- R is a list of repetitions of the uc in the cartesian directions
            - Call self.nd_sup(Q) to build nondiagonal supercell
                -- Q contains the fractional coordinates of the q-point to be folded at Gamma in a nondiagonal supercell like [[m1,m2,m3],[n1,n2,n3]]
            - Call self.displace(modes_file,new_atoms,Temp=0.1) to displace supercell
                -- modes_file is a QE-DFPT phonon-mode output
                -- Temp is the width of the displacements in bohr
                -- new_atoms is the output of (n)d_sup
                
        Sample input script found at:
            tutorials/supercell
            
        NOTA BENE:
            The Q-vector to be folded must be given in CRYSTAL COORDINATES. If computing its related phonons,
            remember that the output of quantum espresso uses CARTESIAN COORDINATES
    """

    def __init__(self,qe_input):

        self.qe_input = qe_input
        self.latvec   = np.array(qe_input.cell_parameters)
        self.basis    = int(qe_input.system['nat'])
        self.atoms    = qe_input.atoms
        self.uc_kpts  = qe_input.kpoints
        self.atypes   = qe_input.atypes
        self.aunits   = qe_input.atomic_pos_type
 
    """
    [START] Displacement-related functions                          
    """
    def displace(self,modes_file,new_atoms,Temp=0.1,use_temp='no',write=True):
        """
        Case of displaced supercell
        """
        #Check if we are displacing the unit cell (i.e., gamma modes)
        GAMMA = False
        try: self.Q
        except AttributeError: GAMMA = True

        print('Applying displacements according to phonon modes...')
        self.use_temp = use_temp
        self.initialize_phonons(modes_file,self.atypes,Temp)
        
        if GAMMA: #No phases and take only optical modes
            phases = np.ones(self.sup_size)
            expand_eigs = np.array([phases[i]*self.eigs for i in range(self.sup_size)])
            self.print_expanded_eigs(expand_eigs,modes_file,GAMMA=GAMMA)
        else: 
            phases = self.getPhases()                
            expand_eigs = np.array([phases[i]*self.eigs for i in range(self.sup_size)])
            self.print_expanded_eigs(expand_eigs,modes_file,GAMMA=GAMMA) #Print expanded eigs
            #Take real part
            for cell in range(self.sup_size): expand_eigs[cell]= self.take_real(expand_eigs[cell])            

        disps = expand_eigs.real.astype(float)
        #Force same gauge choice
        #for cell in range(self.sup_size): disps[cell]= self.force_gauge(disps[cell])
        #Transform eigenmodes in displacements
        #disps[cell][mode][basis][direction]
        disps = np.array([self.osc_length(disp_slice,GAMMA=GAMMA) for disp_slice in disps])
        #disps[mode][cell][basis][direction]
        self.disps = disps.swapaxes(0,1)
        if GAMMA: self.disps = self.disps[3:]    
        if write:
            #A list of PwIn() objects (one for each phonon mode) that can be printed, written to file, etc.
            if 'Q' in globals(): mode='nd'
            else: mode='diagonal'
            self.modes_qe = [self.write(new_atoms,mode,phonon=disps_slice) for disps_slice in self.disps]
    
    def print_expanded_eigs(self,exp_eigs,modes_file,GAMMA=False):
        """
        Print expanded modes in QE-style, to compare and for reference
        """
        if GAMMA:
            q = np.array([0.,0.,0.])
            ncells = 1
            mode_start = 3
        else:
            q = np.array([float(self.Q[0,i])/float(self.Q[1,i]) for i in range(3)])
            ncells = int(np.max(np.array(self.Q[1])))
            mode_start = 0
        
        sc_basis = self.basis*ncells
        freq_THz,freq_cmm1=read_frequencies(modes_file,units='Tera'),read_frequencies(modes_file,units='cmm1')
        exp_eigs = exp_eigs.swapaxes(0,1) #[mode][cell][basis][direction]
        if GAMMA: filename = self.qe_input.control['prefix'][1:-1]+"_s.modes_GAMMA"
        else: filename = self.qe_input.control['prefix'][1:-1]+"_s.modes_expanded"
        exp_eigs_file = open(filename,"w")
        exp_eigs_file.write("     diagonalizing the dynamical matrix ...\n")
        exp_eigs_file.write("\n")
        exp_eigs_file.write(" q =     %2.5f   %2.5f   %2.5f\n"%(q[0],q[1],q[2]))
        exp_eigs_file.write(" **************************************************************************\n") 
        for m in range(mode_start,3*self.basis):
            exp_eigs_file.write("     freq (   %d) =       %2.6f [THz] =       %2.6f [cm-1]\n" \
                                %(m+1,freq_THz[m],freq_cmm1[m]))
            for c in range(ncells):
                for a in range(self.basis):
                    exp_eigs_file.write("( %2.6f  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f )\n" \
                                        %(exp_eigs[m,c,a,0].real,exp_eigs[m,c,a,0].imag, \
                                          exp_eigs[m,c,a,1].real,exp_eigs[m,c,a,1].imag, \
                                          exp_eigs[m,c,a,2].real,exp_eigs[m,c,a,2].imag))
        exp_eigs_file.write(" **************************************************************************")
        exp_eigs_file.close()

    def initialize_phonons(self,modes_file,atypes,Temp):
        #Read frequencies
        Omega = read_frequencies(modes_file)
        self.Omega = Omega*Tera #in hertz
        #Read phonon eigenmodes
        self.eigs = read_eig(modes_file,self.basis)
        #Temperature/displacement
        self.Temp = Temp
        #Atomic masses (in u)
        self.m_at = np.array( [ float( atypes.get( self.atoms[i][0] )[0] ) for i in range(self.basis) ] )

    def getPhases(self):
        q = np.array([float(self.Q[0,i])/float(self.Q[1,i]) for i in range(3)])
        arg = q[0]*self.T[:,0]+q[1]*self.T[:,1]+q[2]*self.T[:,2]
        return np.exp(1j*2.*np.pi*arg)

    def take_real(self,eig):
        modes=3*self.basis
        for i in range(modes):
            #Check that there aren't purely imaginary modes
            if np.real(eig[i]).all == 0.:
                print("mode %i: Purely imaginary"%(i+1))
                eig[i]=1j*eig[i]
        eig = np.real(eig)
        return eig

    def force_gauge(self,eig):
        """ 
        for each normal mode, the first nonzero element is set to be positive
        """
        modes=3*self.basis
        for i in range(modes):
            a=np.nonzero(eig[i])                 # for each normal mode
            a=np.array(a)                        # I get the nonzero elements, of which the first
            if eig[i,a[0,0],a[1,0]] < 0.:        #one, given by these indices, must be >0
                eig[i]=-eig[i]
        return eig

    def osc_length(self,eig,GAMMA=False):
        """
        Oscillator lengths per mode (in ANGSTROM)
        NB:  Only ARBITRARY displacements can be set.
        NB2: Eigenmodes from quantum espresso are ALREADY weighted by atomic masses
        """
        RESCALE     = b2a*self.Temp  #Arbitrary displacement
        temperature = 0.*self.Temp      #Harmonic displacement (permanently set to zero for now)
        modes=3*self.basis
        displacements=[]
        lengths_per_mode=[]
        for nu,eig_slice in enumerate(eig):
            
            if GAMMA and ( nu==0 or nu==1 or nu==2): l=0.
            else: l=sqrt(hbar/(2.*cMp*self.Omega[nu]))
            
            if temperature==0.: sigma2=l*l
            else: sigma2=l*l*(2./(np.exp(hbar*self.Omega[nu]/(kb*temperature))-1.)+1)
            
            if self.use_temp=='no': displacements.append(RESCALE*eig_slice)     #Each mode (i.e. atomic displacement directions) is multiplied by the corresponding length
            else:                   displacements.append(np.sqrt(sigma2)*eig_slice) #Displacement by harmonic sigma
            #List of average "realistic" atomic displacements (not weighted by atomic mass)
            lengths_per_mode.append(np.sqrt(sigma2))
        displacements = np.array(displacements)
        lengths_per_mode = np.array(lengths_per_mode)
        #mass_ratio = np.array([sqrt(Mp/mass) for mass in self.m_at])
        #for d,i in product(range(len(displacements)),range(self.basis)): displacements[d,i,:] *= mass_ratio[i] #Weigh the displ. with the different masses
        
        return displacements #displacements[mode in order of ascending frequency][basis]

    """
    [END] Displacement-related functions   
    """
    
    def lattice_constants(self,vec):
        return [np.linalg.norm(vec[0]),np.linalg.norm(vec[1]),np.linalg.norm(vec[2])]

    def d_sup(self,R,write=True):
        """
        Case of diagonal supercell
        """     
        self.R = R
        self.sup_size = int(R[0]*R[1]*R[2])
        self.new_latvec = np.array([self.latvec[i]*R[i] for i in range(3)])
        self.sup_size = int(self.R[0]*self.R[1]*self.R[2])
        new_atoms = self.build_supercell()
        if write: 
            #PwIn() object that can be printed, written to file, etc.
            self.qe_d = self.write(new_atoms,mode='diagonal')
        return new_atoms

    def nd_sup(self,Q,write=True):
        """
        Case of nondiagonal supercell
        """
        self.Q = np.array(Q)
        print('Nondiagonal supercell')
        if (self.uc_kpts % self.Q[1] != 0).any():
            print('ERROR: You must set a unit cell k-point mesh where%s Nx,Ny,Nz are multiples of %d,%d,%d, respectively.'%('\n',self.Q[1,0],self.Q[1,1],self.Q[1,2]))
            exit()
        self.R, self.new_latvec = self.find_nondiagonal()
        self.sup_size = int(self.R[0]*self.R[1]*self.R[2])
        new_atoms = self.build_supercell()
        if write:
            self.qe_nd = self.write(new_atoms,mode='nd')
        return new_atoms

    def build_supercell(self):
        latvec     = self.latvec
        R          = self.R
        atoms      = np.array([atom[1] for atom in self.atoms])
        if self.aunits!='angstrom': atoms = red_car(atoms,latvec) 
        else: latvec = b2a*latvec
        #new_atoms[cell][basis][direction]
        new_atoms      = np.array([atoms for n in range(self.sup_size)])
        T = []
        for nz,ny,nx in product(range(int(R[2])),range(int(R[1])),range(int(R[0])) ): 
            cell=int(nx+ny*R[0]+nz*R[0]*R[1]) 
            translation = nx*latvec[0] +ny*latvec[1] +nz*latvec[2]
            for b in range(self.basis): 
                new_atoms[cell,b]=new_atoms[cell,b] + translation 
            T.append(translation)
        T = np.array(T) #Positions of the repeated unit cells
        if self.aunits!='angstrom': self.T=car_red(T,self.latvec)
        else: self.T=T
        #new_atoms[super_basis][directions]$
        new_atoms=new_atoms.reshape(self.basis*self.sup_size,3) 
        if self.aunits!='angstrom': new_atoms = car_red(new_atoms,self.new_latvec)
        return new_atoms

    def find_integers(self,nums,g23,g12,g31,g123):
        """
        Compute integers for off-diagonal supercell matrix elements 
        Called by find_nondiagonal()
        """
        if nums[1]==0: p=0
        else: 
            #Compute p (it's a modulo equation)
            if g23 == 1: p = 0
            else:
                for i in range(1,g23):
                    if (nums[1]+i*nums[2]) % g23 == 0:
                        p=i
                        break
        if nums[0]==0: q,r=[0,1] #[POSSIBLE BUG for certain q-vectors] These conditions must be checked carefully
        else:
            #Compute q
            g12_r = int(g12/g123)
            g23_r = int(g23/g123)
            g31_r = int(g31/g123)
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
        """
        Nondiagonal supercell, based on [Phys. Rev. B 92, 184301]
        """
        Q = self.Q
        #Take care of components already at Gamma
        Q[1,np.where(Q[0]==0)]=1
        #Shift the q-point into the positive quadrant of the reciprocal unit cell
        Q[0,np.where(Q[0]<0)]+=Q[1,np.where(Q[0]<0)]
        #GCDs of Q[1] (in the logical order of the derivation)
        g23  = gcd(Q[1,1],Q[1,2])
        g12  = gcd(Q[1,0],Q[1,1])
        g31  = gcd(Q[1,2],Q[1,0])
        g123 = gcd(Q[1,0],gcd(Q[1,1],Q[1,2]))
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
        """
        Function to compute reciprocal lattice
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
        """ 
        Put the atomic element labels in the right order
        """
        positions_input = new_atoms.tolist()
        elements_input  = [[self.qe_input.atoms[i][0] for i in range(int(self.basis))] for j in range(self.sup_size)]  
        elements_input  = [ item for sublist in elements_input for item in sublist ]
        atoms_input     = [[elements_input[i], positions_input[i]] for i in range(self.sup_size*self.basis)]
        return atoms_input

    def posint(self,value):
        return abs(int(round(value)))

    def write(self,new_atoms,mode,phonon=None):
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
        if phonon is not None:
            phonon = phonon.reshape(self.basis*self.sup_size,3)
            phonon = car_red(phonon,self.new_latvec)
            new_atoms = new_atoms + phonon
        qe_s = copy.deepcopy(qe)
        qe_s.set_atoms(self.atoms_input(new_atoms))
        qe_s.control['prefix'] = qe.control['prefix'][:-1]+"_s'"
        #[POSSIBLE BUG] with only ibrav==0 and cell_parameters, it might fail the symmetry !?
        if 'celldm(1)' in qe_s.system: del qe_s.system['celldm(1)']
        if 'celldm(2)' in qe_s.system: del qe_s.system['celldm(2)']
        if 'celldm(3)' in qe_s.system: del qe_s.system['celldm(3)']
        """
        qe_s.system['celldm(1)'] = alat[0]
        qe_s.system['celldm(2)'] = alat[1]/alat[0]
        qe_s.system['celldm(3)'] = alat[2]/alat[0]
        """
        qe_s.system['ibrav']=0
        qe_s.cell_units = 'bohr'
        qe_s.cell_parameters = new_latvec
        #Just a suggestion for the new bands
        if 'nbnd' in qe.system: qe_s.system['nbnd'] = self.sup_size*int(qe.system['nbnd'])
        qe_s.system['nat'] = self.basis*self.sup_size
        qe_s.kpoints = new_kpoints
        return qe_s


