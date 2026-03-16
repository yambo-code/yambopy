import numpy as np
import sys
from qepy import *

from yambopy.units import *
from yambopy.zeros import *
from yambopy.tools.funcs import bose

import os

# convertion constants
eVtocm1 = 8065.54429
cm1toeV = 1.0/eVtocm1
ha2ev  = 27.211396132
eV2ha  = 1.0/ha2ev
Thz2cm1 = 33.35641
cm12Thz = 1.0/33.35641

"""
Starting from a unit cell (uc) input that is read from file, we build:

- diagonal supercell
- non-diagonal supercell folding a specific q-point Q
- displaced non-diagonal supercells along phonon eigenmodes at Q (read from file)
"""

def generate_diagonal_supercell(uc,R):
    """
    First case: diagonal supercell of size R
    """
    qe = PwIn.from_file(uc) #read the uc input

    sc = Supercell(qe) # initialize class
    sc.d_sup(R) # Generate supercell as PwIn object called 'qe_d'

    #name of output file
    suffix = ''.join(sc.qe_d.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_d.control['prefix'].split('\''))

    #write supercell to file
    sc.qe_d.write('sc_diagonal_%s.%s'%(prefix,suffix))

    print('Diagonal supercell written to file.')

def generate_nondiagonal_supercell(uc,Q,kpoints=None):
    """
    Second case: nondiagonal supercell folding point Q
    """
    qe = PwIn.from_file(uc) #read the uc input
    if kpoints is not None: qe.kpoints = kpoints #Optional: manually change the original kpt mesh to ensure consistency    

    sc = Supercell(qe) # initialize class
    sc.nd_sup(Q) # Generate supercell as PwIn object called 'qe_nd'

    #name of output file
    suffix = ''.join(sc.qe_nd.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_nd.control['prefix'].split('\''))

    #write supercell to file
    sc.qe_nd.write('sc_nondiagonal_%s.%s'%(prefix,suffix))

    print('Nondiagonal supercell written to file.')   

def generate_displaced_supercells(uc,Q,modes_file,kpoints=None, use_temp = 'no', temp = 0.1):
    """
    Third case: displaced supercells along phonon modes at Q
    """
    qe = PwIn.from_file(uc) #read the uc input
    if kpoints is not None: qe.kpoints = kpoints #Optional: manually change the original kpt mesh to ensure consistency

    sc = Supercell(qe) # initialize class
    nd_atom_positions = sc.nd_sup(Q) # Generate supercell as PwIn object called 'qe_nd' getting new atomic positions
    
    #Read the matdyn file
    qe_dyn=Matdyn.from_modes_file(filename=modes_file)

    # Displace atoms./
    # Intensity is Temp (in bohr)
    # Sign and direction (standing wave at Q) given by modes_file
    #sc.displace(modes_file,nd_atom_positions,Temp=0.1) # Generate list of displaced supercells as PwIn objects called 'modes_qe'
    # iq index of the q-point in the matdyn file (default 0, the first q-point)
    sc.displace(qe_dyn,nd_atom_positions,iq=0,Temp=temp, use_temp=use_temp)
    N_modes = len(sc.modes_qe)

    #name of output file
    suffix = ''.join(sc.qe_nd.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_nd.control['prefix'].split('\''))

    #write supercells to file for each mode
    for mode in range(N_modes): sc.modes_qe[mode].write('sc_displaced_%s_mode%d.%s'%(prefix,mode+1,suffix)) 

    print('Displaced supercells written to file.')

def generate_displaced_unitcell(uc,modes_file):
    """
    Fourth case: displaced cell at Q=0
    """
    qe = PwIn.from_file(uc) #read the uc input

    sc = Supercell(qe) # initialize class
    atom_positions = sc.d_sup([1,1,1]) # Generate "supercell" with size 1 as PwIn object called 'qe_d' getting atomic positions
    
    #Read the matdyn file
    qe_dyn=Matdyn.from_modes_file(filename=modes_file)

    # Displace atoms.
    # Intensity is Temp (in bohr)
    # Sign and direction (standing wave at Q) given by modes_file
    sc.displace(qe_dyn,atom_positions,Temp=0.1) # Generate list of displaced supercells as PwIn objects called 'modes_qe'
    N_modes = len(sc.modes_qe)

    #name of output file
    suffix = ''.join(sc.qe_d.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_d.control['prefix'].split('\''))

    #write supercells to file for each mode
    for mode in range(N_modes): sc.modes_qe[mode].write('uc_GAMMA_displaced_%s_mode%d.%s'%(prefix,mode+4,suffix))

def sort_all_phonon_modes(matdyn:Matdyn):
    print('Sort phonon modes of whole ensemble, excluding COM motion modes...') 
    eigvals = matdyn.eig.flatten() # indexed by i = i(iq, imode) = iq * nmodes + imode
    arg_sort = np.argsort(eigvals)
    ntot_modes = matdyn.nqpoints * matdyn.nmodes
    # eigvecs = np.zeros([ntot_modes, matdyn.natoms * 3], dtype=complex)
    eigvecs = np.reshape(matdyn.eiv, [ntot_modes, matdyn.natoms * 3])
    qlist = np.zeros([ntot_modes, 3], dtype=float)
    i = 0
    for iq, imode in product(np.arange(matdyn.nqpoints), np.arange(matdyn.nmodes)):
        # i = iq * matdyn.nmodes + imode
        # eigvecs[i] = matdyn.eiv[iq, imode, :]
        qlist[i] = matdyn.qpoints[iq]
        i+= 1
    eigvals = eigvals[arg_sort][3:]
    eigvecs = eigvecs[arg_sort][3:]
    qlist   = qlist[arg_sort][3:]
    return qlist, eigvals, eigvecs

def generate_thermal_displaced_supercells(uc, Q, modes_file, kpoints=None, T = 0., qe_control_dict=None, qe_system_dict=None):
    """
    Thermal displacement at a fixed Q of phonon. 
    """
    qe = PwIn.from_file(uc) #read the uc input
    if kpoints is not None: qe.kpoints = kpoints #Optional: manually change the original kpt mesh to ensure consistency
    sc = MySupercell(qe) # initialize class
    nd_atom_positions = sc.nd_sup(Q) # Generate supercell as PwIn object called 'qe_nd' getting new atomic positions
    #Read the matdyn file
    qe_dyn=Matdyn.from_modes_file(filename=modes_file)
    # displace atoms in each supercell according to the mode and temperature
    sc.displace(qe_dyn, nd_atom_positions, iq=0, T = T)
    N_modes = len(sc.modes_qe)
    #name of output file
    suffix = ''.join(sc.qe_nd.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_nd.control['prefix'].split('\''))
    #write supercells to file for each mode
    for mode in range(N_modes): 
        if qe_control_dict is not None:
            for k, v in qe_control_dict.items():
                sc.modes_qe[mode].control[k] = v
        if qe_system_dict is not None:
            for k, v in qe_system_dict.items():
                sc.modes_qe[mode].system[k] = v
        sc.modes_qe[mode].write('sc_displaced_%s_mode%d.%s'%(prefix,mode+1,suffix)) 
    print('Displaced supercells written to file.')


class MySupercell(Supercell):
    def __init__(self, qe_input:PwIn):
        super().__init__(qe_input)

    def initialize_phonons(self,iq:int, qe_dyn:Matdyn):
        '''
        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Overloads the original function with Temp parameter removed. 
        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        '''
        self.qe_dyn   = qe_dyn
        self.iq       = iq
        self.Omega    = qe_dyn.eig[iq]*Tera
        self.eiv = np.reshape(self.qe_dyn.eiv[iq],(self.qe_dyn.nmodes,self.basis,3))

    def displace(self, qe_dyn:Matdyn, new_atoms, iq=0, T = 0.0, uniform_disp = 0.1, use_thermal_disp = True, write=True): # T is explicitly for Temperature, in Kelvin.
        """
        Case of displaced supercell
        """
        #Check if we are displacing the unit cell (i.e., gamma modes)
        GAMMA = False
        try: self.Q
        except AttributeError: GAMMA = True

        if use_thermal_disp: 
            print('Thermal displacement activated. Will ignore *uniform_disp*.')
        
        print('Applying displacements according to phonon modes...')
        self.initialize_phonons(iq, qe_dyn)

        if GAMMA and np.linalg.norm(qe_dyn.qpoints[iq])>matdyn_q_atol:
            print("WARNING: Q-point in matdyn file different from the one used to generate supercell")
        if not GAMMA:
            # Transform all q-point in reduced coordinates
            rlat =rec_lat(self.latvec)
            alat0=self.qe_input.get_alat0()
            q = []
            q.append(qe_dyn.qpoints[iq]/alat0)
            q_matdyn=car_red(q,rlat)
            q_in = np.array([float(self.Q[0,i])/float(self.Q[1,i]) for i in range(3)])

            # Bring in the BZ [0,1)
            q_in = q_in - np.floor(q_in)
            q_matdyn = q_matdyn - np.floor(q_matdyn)
            if not np.allclose(q_matdyn,q_in,rtol=0.0,atol=matdyn_q_atol):
                print("WARNING: Q-point in matdyn file different from the one used to generate supercell")
                print("Q-in     : ",q_in, " [red] ")
                print("Q-matdyn : ",q_matdyn[0], " [red] ")

        if GAMMA: #No phases and take only optical modes
            phases = np.ones(self.sup_size)
            expand_eigs = np.array([phases[i]*self.eiv for i in range(self.sup_size)])
            self.print_expanded_eigs(expand_eigs,GAMMA=GAMMA)
        else: 
            phases = self.getPhases()                
            expand_eigs = np.array([phases[i]*self.eiv for i in range(self.sup_size)])
            self.print_expanded_eigs(expand_eigs,GAMMA=GAMMA) #Print expanded eigs
            #Take real part
            for cell in range(self.sup_size): expand_eigs[cell]= self.take_real(expand_eigs[cell])            

        disps = expand_eigs.real.astype(float) # disps here are actually the eigen vectors (real parts). 
        #Force same gauge choice
        #for cell in range(self.sup_size): disps[cell]= self.force_gauge(disps[cell])
        #Transform eigenmodes in displacements
        #disps[cell][mode][basis][direction]
        disps = np.array([self.osc_length(disp_slice, T, uniform_disp, use_thermal_disp, GAMMA=GAMMA) for disp_slice in disps])
        #disps[mode][cell][basis][direction]
        self.disps = disps.swapaxes(0,1)
        if GAMMA: self.disps = self.disps[3:]    
        if write:
            #A list of PwIn() objects (one for each phonon mode) that can be printed, written to file, etc.
            if 'Q' in globals(): 
                mode='nd'
            else: 
                mode='diagonal'
            self.modes_qe = [self.write(new_atoms,mode,phonon=disps_slice) for disps_slice in self.disps]
    
    def osc_length(self, eig, T = 0.0, uniform_disp = 0.1, use_thermal_disp = True, GAMMA=False): # Enables temperature
        """
        +++++++++++++++++++++++++++++++++++++++++++
        Overloads the original osc_length() method.
        +++++++++++++++++++++++++++++++++++++++++++
        Oscillator lengths per mode (in ANGSTROM)
        NB: If *use_thermal_disp* is set to True, uniform_disp will be ignored. -- YUNCHENG MAO
        NB2: Eigenmodes from quantum espresso are ALREADY weighted by atomic masses.

        TODO:
        Displacement by Gamma phonons not yet well done. Need to fix this. 

        MEM: Why GAMMA should be treated separately??? I don't see the point...
        One specifies the supercell, and calculates the phases according to phonon Q. GAMMA should be treated the same.
        """
        # RESCALE     = b2a * self.Temp  #Arbitrary displacement
        # temperature = 0.0 * self.Temp      #Harmonic displacement (permanently set to zero for now)
        # modes=3*self.basis
        displacements=[]
        # lengths_per_mode=[]
        for nu, eig_slice in enumerate(eig):
            if use_thermal_disp:
                if GAMMA and ( nu==0 or nu==1 or nu==2): 
                    '''
                    It tries to throw aways the 3 lowest energy phonons at Gamma, the three collective modes of the whole cell.
                    i.e. the movement is all in the center of mass.
                    '''
                    q_T=0.0
                # w_au = self.qe_dyn.get_phonon_freq(0, nu+1, unit='Ha') # Stupidly written function. 
                w_au = self.qe_dyn.eig[0, nu]*cm1toeV*eV2ha # NB: having only single phonon Q.
                q_0  = 1.0/math.sqrt(2.0 * w_au) # Is this factor 2 correct ????
                q_T  = q_0*math.sqrt(1.0 + 2.0 * bose(w_au * ha2ev, T))
                print(f"Zero-temperature amplitude for mode {nu+1} is {q_0:.8g} (a.u.).")
                print(f"At given Temperature ({T:.2g} K), amplitude for mode {nu+1} is {q_T:.8g} (a.u.).")
                displacements.append(q_T * eig_slice/np.sqrt(amu2au))
            else:
                displacements.append(b2a * uniform_disp * eig_slice)
            # else: 
            #     l=sqrt(hbar/(2.*cMp*self.Omega[nu]))
            # if temperature==0.: 
            #     sigma2=l*l
            # else: 
            #     sigma2=l*l*(2./(np.exp(hbar*self.Omega[nu]/(kb*temperature))-1.)+1)
            # if self.use_temp=='no': 
            #     print(f"Manual uniform displacement by {self.Temp:.4g} bohr along each eigen mode.")
            #     displacements.append(b2a * self.Temp * eig_slice)     #Each mode (i.e. atomic displacement directions) is multiplied by the corresponding length
            # else:                   
            #     displacements.append(q_T * eig_slice) 
            #List of average "realistic" atomic displacements (not weighted by atomic mass)
            # lengths_per_mode.append(np.sqrt(sigma2))
        displacements = np.array(displacements)
        # lengths_per_mode = np.array(lengths_per_mode)
        #mass_ratio = np.array([sqrt(Mp/mass) for mass in self.m_at])
        #for d,i in product(range(len(displacements)),range(self.basis)): displacements[d,i,:] *= mass_ratio[i] #Weigh the displ. with the different masses
        
        return displacements #displacements[mode in order of ascending frequency][basis]


    def apply_full_thermal_displacement(self, R_sc:np.ndarray, matdyn:Matdyn, T_Kelvin:float, 
                                        rm_output_sym = False,write_wf=True,force_symmorphic=True,
                                        qe_control_dict:dict = {}, qe_system_dict:dict = {}, 
                                        qe_electron_dict:dict={}, output_all=False, temp_dir = None):
        '''
        R_sc: e.g. an integer array [3, 3, 1] defining the DIAGONAL supercell 
        matdyn: the full phonon output from normal unit cell.
        ''' 
        '''
        Unormalize vector if necessary 
        '''
        if not matdyn.check_orthogonality():
            print("Eigenvectors not normlized unscale masses")
            matdyn.unnormalize_with_masses(self.qe_input.get_masses())
        else:
            print("Eigenvectors are normlized")

        qlist, eigvals, eigvecs = sort_all_phonon_modes(matdyn) # Now every thing is ordered in ascending order of phonon energy.
        '''
        Note: phonon q in qlist are in units of 2 pi/ alat
        '''
        ntot_modes = len(eigvals)
        print(f'INFO: {ntot_modes} modes in total.')
        #
        # Check if number of q-points is compatible with cell size
        #
        ntot_cells=np.prod(R_sc)
        if ntot_cells != matdyn.nqpoints:
            raise ValueError('Number of q-points not compatible with supercell ')
        #
        # Check that q-vectors are compatible with the cell
        #
        q = []
        for nq in range(matdyn.nqpoints):
            q.append(matdyn.qpoints[nq]/self.qe_input.get_alat0()) # in reduced units of reciprocal lattice bases
        q_red=car_red(q,rec_lat(self.latvec))
        for nq in range(matdyn.nqpoints):
            q_red[nq]=q_red[nq,:]*R_sc[:]
        if not np.all(np.isclose(q_red, np.round(q_red), rtol=0, atol=matdyn_q_atol)):
        #    print(q_red)
            raise ValueError('Q-vectors not compatible with the supercell')
        else:
            print("Q-vectors compatible with the supercell :",R_sc)

        self.full_eigvecs = eigvecs.reshape([ntot_modes, matdyn.natoms, 3])
        eigvals = eigvals * cm1toeV * eV2ha # convert phonon energy to Ha
        # displacement amplitudes are only determined by phonon energy and temperature
        disp_amplitudes = 1.0/np.sqrt(2.0 * eigvals)
        disp_amplitudes *= np.sqrt(1.0 + 2.0 * bose(eigvals * ha2ev, T_Kelvin))/np.sqrt(amu2au)
        self.sc_atom_positions = self.d_sup(R_sc) # return non-displaced atom coordinates in the DIAGONAL supercell
        newlat_constants = self.lattice_constants(self.new_latvec)
        self.qe_d.system['celldm(1)'] = newlat_constants[0]
        self.qe_d.system['celldm(2)'] = newlat_constants[1]/newlat_constants[0]
        self.qe_d.system['celldm(3)'] = newlat_constants[2]/newlat_constants[0]
        Qs = car_red(qlist, rec_lat(np.array(self.qe_input.cell_parameters)/self.qe_input.get_alat0())) # in reduced units of reciprocal lattice bases
        # one line of code to replace the getPhases() function, and does more things.
        phases = np.exp(2j*np.pi * self.T @ Qs.T)# resulted phases are indexed as: phase[icell, imode] 
        phases = phases.T # now indexed as phase[imode, icell]
        '''
        Next: expand eigenvectors, and compute displacements. 
        NB: AVOID using loops as much as possible!!! 
        Carefully align the matrix elements and let the low-level code do the math.
        '''
        self.expanded_phases = [] # for debug
        # expand_eigvecs = np.zeros([ntot_modes, self.sup_size * matdyn.natoms * 3], dtype=complex) # Eigenvectors represented in coordiantes of ALL ATOMS IN THE SUPERCELL
        expanded_phases = np.zeros([ntot_modes, self.sup_size * self.basis * 3], dtype=complex)
        expanded_eigvecs = np.hstack([eigvecs] * self.sup_size)
        for i, ph_mode in enumerate(phases):
            expanded_phases[i] = np.hstack([ph_mode.reshape([-1,1])] * (self.basis * 3)).flatten()
            # phases_to_apply = np.hstack([phases[i].reshape(-1,1)] * (matdyn.natoms * 3)).flatten() # repeat phases to every coordinate in each replicated cell
            # expand_eigvecs[i] = np.hstack([eigvecs[i]] * self.sup_size) * phases_to_apply
            # self.expanded_phases.append(phases_to_apply)
        expanded_eigvecs *= expanded_phases
        '''
        Normalize eigenvectors
        '''
        for i in range(len(expanded_eigvecs[:,0])):
            expanded_eigvecs[i]/=np.linalg.norm(expanded_eigvecs[i])
        '''
        Rescale by masses
        '''
        qe_non_disp = self.write(self.sc_atom_positions, mode='diagonal')
        masses_s=qe_non_disp.get_masses()
        for n in range(len(expanded_eigvecs[:,0])):
            for a in range(qe_non_disp.natoms):
                        expanded_eigvecs[i,a*3:(a+1)*3] *= 1.0/sqrt(masses_s[a])

        # A list of displacements
        # expand_disp_amplitudes = np.hstack([disp_amplitudes.reshape(-1,1)] * (self.sup_size * matdyn.natoms * 3))
        # displacements = np.real(expand_disp_amplitudes * expanded_eigvecs) # displacements[imode, ix]
        displacements = np.real(np.diag(disp_amplitudes) @ expanded_eigvecs)
        if output_all:
            if temp_dir is None: 
                    temp_dir = '_all_inputs'
            if not os.path.exists(temp_dir):
                os.mkdir(temp_dir)
            np.savetxt(os.path.join(temp_dir, 'Q_coors.txt'), Qs, fmt='%7.4f')
            np.savetxt(os.path.join(temp_dir, 'Q_carts.txt'), qlist, fmt='%7.4f')
            for i, disp in enumerate(displacements):
                qe_s = self.write(self.sc_atom_positions, 'diagonal', disp)
                qe_s.write(os.path.join(temp_dir, f'mode_{i+1}.in'))
                
        # tot_disp = np.sum(np.diag((-1)**np.arange(ntot_modes)) @ displacements, axis=0)
#        tot_disp = np.sum(displacements, axis=0)/matdyn.nqpoints
        tot_disp = np.sum(displacements, axis=0)/matdyn.nqpoints
        
        '''
        NB:
        DO NOT use the displace() method from PwIn class. The mass weights are already applied to the phonon modes in the QE output. 
        This function wrongly multiplies once again the mass weights.
        '''
        # Take advantage of Supercell.write() to assign the displacements to the supercell atoms
        # Then use PwIn.write() to save the QE input files.
        for i in range(2):
            qe_s = self.write(self.sc_atom_positions, 'diagonal', (-1)**i * tot_disp)
            qe_s.control['pseudo_dir'] = "'./'"
            if force_symmorphic:
                qe_s.system['force_symmorphic'] = '.true.' # important for later Yambo calculations
            for k, v in qe_control_dict.items():
                qe_s.control[k] = v
            for k, v in qe_system_dict.items():
                qe_s.system[k] = v
            #---- output to files ----
            qe_s.control['calculation'] = "'scf'"
            fname = 'thermal_disp_%d.scf.in' % (i+1)
            print('Writing file %s ...' % fname)
            qe_s.write(fname)
            qe_s.control['calculation'] = "'nscf'"
            if rm_output_sym:
                qe_s.system['nosym'] = '.true.'
                qe_s.system['noinv'] = '.true.'
            if not write_wf:
                qe_s.control['disk_io']= "'none'"
            if 'nbnd' not in qe_s.system.keys():
                print('Warning: nbnd not specified!')
            fname = 'thermal_disp_%d.nscf.in' % (i+1)
            print('Write file %s ...' % fname)
            qe_s.write(fname)
        #----------------- debug purpose -----------------------
        self.phases = phases
        self.expanded_eigvecs = expanded_eigvecs
        self.expanded_phases = np.array(self.expanded_phases)
        self.Qs = Qs
        self.displacements = displacements
        self.disp_amplitudes = disp_amplitudes
        self.tot_disp = tot_disp
        #----------------- debug purpose -----------------------
    
    def gen_displaced_supercells(self, iq:int, matdyn:Matdyn, mode:str, T_Kelvin:float = 0, sc_spec_dict:dict = {}, qe_control_dict:dict={}, qe_system_dict:dict={}, rm_sym=True):
        # generate a supercell according to a single q of phonon. 
        q = matdyn.qpoints[iq]
        if mode in ['diagonal', 'd']:
            mode = 'diagonal' # to be compatible with parameters in Supercell.write()
            new_atom_positons = self.d_sup(sc_spec_dict['R'])
            # Q = 1/np.array(sc_spec_dict['R'])
            Q = car_red([q], rec_lat(np.array(self.qe_input.cell_parameters)/self.qe_input.get_alat0())).flatten()
            newlat_constants = self.lattice_constants(self.new_latvec)
            self.qe_d.system['celldm(1)'] = newlat_constants[0]
            self.qe_d.system['celldm(2)'] = newlat_constants[1]/newlat_constants[0]
            self.qe_d.system['celldm(3)'] = newlat_constants[2]/newlat_constants[0]
        elif mode in ['nondiagonal', 'nd']:
            mode = 'nondiagonal' # to be compatible with parameters in Supercell.write()
            new_atom_positons = self.nd_sup(sc_spec_dict['Q'])
            q_red = car_red([q], rec_lat(np.array(self.qe_input.cell_parameters)/self.qe_input.get_alat0())).flatten()
            Q = np.array(sc_spec_dict['Q'])
            Q = Q[0]/Q[1]
            if np.linalg.norm(q_red - Q) > 0.1:
                print("Warning: phonon q specified by sc_spec_dict['Q'] is very different from the one selected with iq.")
            newlat_constants = self.lattice_constants(self.new_latvec)
            self.qe_nd.system['celldm(1)'] = newlat_constants[0]
            self.qe_nd.system['celldm(2)'] = newlat_constants[1]/newlat_constants[0]
            self.qe_nd.system['celldm(3)'] = newlat_constants[2]/newlat_constants[0]
        else:
            print("ERROR: Unrecognized mode specification. Should be one of 'd', 'diagonal', 'nd' and 'nondiagonal'. ")
            return 
        eigvals = matdyn.eig[iq] * cm1toeV * eV2ha
        eigvecs = matdyn.eiv[iq]
        nmodes = matdyn.nmodes
        # Remove COM modes from Gamma phonon
        if np.linalg.norm(Q) < 0.001: 
            print('This is GAMMA phonon. Remove COM modes.')
            eigvals = eigvals[3:]
            eigvecs = eigvecs[3:]
            nmodes = nmodes - 3
        disp_amplitudes = 1/np.sqrt(2 * eigvals)
        disp_amplitudes *= np.sqrt(1.0 + 2.0 * bose(eigvals * ha2ev, T_Kelvin))/np.sqrt(amu2au)
        phases = np.exp(2j* np.pi * self.T @ Q)
        # expanded_eigvecs = np.zeros([nmodes, self.sup_size * matdyn.natoms * 3], dtype=complex)
        expanded_phases = np.hstack([phases.reshape([-1,1])] * (matdyn.natoms * 3)).flatten()
        expanded_eigvecs = np.hstack([eigvecs] * self.sup_size) * expanded_phases
        displacements = np.real(np.diag(disp_amplitudes) @ expanded_eigvecs)
        # displacements_red = car_red(displacements, self.new_latvec)
        for imode in np.arange(nmodes):
            print(f'Writing mode {imode+1}')
            qe_s = self.write(new_atom_positons, mode=mode, phonon=displacements[imode])
            qe_s.control['pseudo_dir'] = "'./'" # fixes the bug from qepy
            qe_s.system['force_symmorphic'] = '.true.'
            for k, v in qe_control_dict.items(): 
                qe_s.control[k] = v
            for k, v in qe_system_dict.items():
                qe_s.system[k] = v
            qe_s.control['calculation'] = "'scf'"
            qe_s.write('mode_%d.scf.in' % (imode+1))
            qe_s.control['calculation'] = "'nscf'"
            if rm_sym:
                qe_s.system['nosym'] = '.true.'
                qe_s.system['noinv'] = '.true.'
            qe_s.write('mode_%d.nscf.in' % (imode + 1))
        # debug
        self.new_atom_positions = new_atom_positons
        self.displacements = displacements

        
