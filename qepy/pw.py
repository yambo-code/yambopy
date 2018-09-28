# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
import re
import shutil
from math import sqrt
from .pseudo import get_pseudo_path
from .tools import fortran_bool

class PwIn(object):
    """
    Class to generate an manipulate Quantum Espresso input files
    Can be initialized either reading from a file or starting from a new file.
    This class is not meant to be comprehensive but a lightweight version capable of basic input/ouput of PW files.
    For a comprehensive class use the tools provoded by the AiiDa package: http://www.aiida.net/

    Examples of use:

    To read a local file with name "mos2.in"

        .. code-block :: python

            qe = PwIn.from_file('mos2.scf')
            print qe

    To start a file from scratch

        .. code-block :: python

            qe = PwIn('mos2.scf')
            qe.atoms = [['N',[ 0.0, 0.0,0.0]],
                        ['B',[1./3,2./3,0.0]]]
            qe.atypes = {'B': [10.811, "B.pbe-mt_fhi.UPF"],
                         'N': [14.0067,"N.pbe-mt_fhi.UPF"]}

            qe.control['prefix'] = "'%s'"%prefix
            qe.control['verbosity'] = "'high'"
            qe.control['wf_collect'] = '.true.'
            qe.control['pseudo_dir'] = "'../pseudos/'"
            qe.system['celldm(1)'] = 4.7
            qe.system['celldm(3)'] = layer_separation/qe.system['celldm(1)']
            qe.system['ecutwfc'] = 60
            qe.system['occupations'] = "'fixed'"
            qe.system['nat'] = 2
            qe.system['ntyp'] = 2
            qe.system['ibrav'] = 4
            qe.kpoints = [6, 6, 1]
            qe.electrons['conv_thr'] = 1e-8

            print qe

    Special care should be taken with string variables e.g. "'high'"

    """
    _pw = 'pw.x'

    def __init__(self):
        """ TODO: specify the required parameters """
        #kpoints
        self.ktype = "automatic"
        self.kpoints = [1,1,1]
        self.shiftk = [0,0,0]
        self.klist = []
        #dictionaries
        self.control = dict(prefix="'pw'",wf_collect='.true.')
        self.system = dict()
        self.electrons = dict(conv_thr=1e-8)
        self.ions = dict()
        self.cell = dict()
        self.atypes = dict()
        self.atoms = []
        self.cell_parameters = []
        self.cell_units = 'angstrom'
        self.atomic_pos_type = 'crystal'

    @classmethod
    def from_file(cls,filename):
        """ Initialize the QE structure from a file """
        new = cls()
 
        with open(filename,"r") as f:
            new.file_lines = f.readlines() #set file lines
            new.store(new.control,"control")     #read &control
            new.store(new.system,"system")      #read &system
            new.store(new.electrons,"electrons")   #read &electrons
            new.store(new.ions,"ions")        #read &ions
            new.store(new.cell,"cell")        #read &ions
            #read ATOMIC_SPECIES
            new.read_atomicspecies()
            #read ATOMIC_POSITIONS
            new.read_atoms()
            #read K_POINTS
            new.read_kpoints()
            #read CELL_PARAMETERS
            new.read_cell_parameters()

        return new

    @classmethod
    def from_structure_dict(cls,structure,kpoints=None,ecut=None,pseudo_dir='.'):
        pwi = cls()
        pwi.set_structure(structure)
        if kpoints: pwi.set_kpoints(kpoints)
        if ecut: pwi.set_ecut(ecut)
        if pseudo_dir: pwi.pseudo_dir = pseudo_dir
        return pwi
       
    @property
    def pseudo_dir(self):
        return self.control['pseudo_dir'].replace("'",'')
    
    @pseudo_dir.setter
    def pseudo_dir(self,value):
        self.control['pseudo_dir'] = "'%s'"%value.replace("'",'')

    @property
    def prefix(self): 
        return self.control['prefix'].replace("'",'')

    @prefix.setter
    def prefix(self,value):
        self.control['prefix'] = "'%s'"%value.replace("'",'')

    def set_ecut(self,ecut):
        self.system['ecutwfc'] = ecut

    def set_structure(self,structure):
        """
        Set the structure from a structure dictionary
        
        Example:
            .. code-block :: python
            
            structure = dict(ibrav=4,celldm1=4.7,celldm3=12)
            atypes = dict(Si=[28.086,"Si.pbe-mt_fhi.UPF"])
            atoms = [['N',[ 0.0, 0.0,0.5]],
                     ['B',[1./3,2./3,0.5]]]
            BN = dict(structure=structure,atypes=atypes,atoms=atoms)

            pwi = PwIn()
            pwi.set_structure(BN)
        """
        if 'lattice' in structure: self.set_lattice(**structure['lattice'])
        if 'atypes'  in structure: self.set_atypes(structure['atypes'])
        if 'atoms'   in structure: self.set_atoms(structure['atoms'])

    def set_lattice(self,ibrav=None,celldm1=None,celldm2=None,celldm3=None,
                      celldm4=None,celldm5=None,celldm6=None):
        """Set the structure using the typical QE input variables"""
        if ibrav is not None: self.system['ibrav'] = ibrav 
        if celldm1 is not None: self.system['celldm(1)'] = celldm1
        if celldm2 is not None: self.system['celldm(2)'] = celldm2
        if celldm3 is not None: self.system['celldm(3)'] = celldm3
        if celldm4 is not None: self.system['celldm(4)'] = celldm4
        if celldm5 is not None: self.system['celldm(5)'] = celldm5
        if celldm6 is not None: self.system['celldm(6)'] = celldm6

    def set_atoms(self,atoms):
        """
        Set the atoms
        
        Example:
            .. code-block :: python
   
                pwi = PwIn()
                pwi.set_atoms( [['Si',[0.125,0.125,0.125]],
                                ['Si',[-.125,-.125,-.125]]])
  
        """
        self.system['nat'] = len(atoms)
        #TODO: add consistency check
        self.atoms = atoms

    def set_atypes(self,atypes):
        """"
        Set the atom types.

        Example:
            .. code-block :: python
                pwi = PwIn()
                pwi.set_atypes({'Si': [28.086,"Si.pbe-mt_fhi.UPF"]})
        """
        self.system['ntyp'] = len(atypes)
        #TODO: add consistency check
        self.atypes = atypes

    def set_nscf(self,nbnd,conv_thr=1e-8,diago_full_acc=True,force_symmorphic=True):
        """
        set the calculation to be nscf
        """
        self.control['calculation'] = "'nscf'"
        self.electrons['conv_thr'] = conv_thr
        self.system['nbnd'] = nbnd
        self.electrons['diago_full_acc'] = fortran_bool(diago_full_acc)
        self.system['force_symmorphic'] = fortran_bool(force_symmorphic)
        return self

    def get_pseudos(self,destpath='.',pseudo_paths=[],verbose=0):
        """
        Check if the pseudopotential files can be found in the specified path
        and copy them to destpath
        """
        import qepy.data.pseudos as qe_pseudos
        pseudo_paths.append(os.path.dirname(qe_pseudos.__file__))

        #use pseudo_dir from control
        if 'pseudo_dir' in self.control:
            pseudo_dir = self.control['pseudo_dir']
            if os.path.isdir(pseudo_dir): pseudo_paths.append(pseudo_dir)
        
        ppstring = '\n'.join(pseudo_paths)
        if verbose: print('List of pseudo_paths:\n'+ppstring)

        #find all the pseudopotentials
        for atype,(mass,pseudo_filename) in self.atypes.items():
            pseudo_filename = get_pseudo_path( pseudo_filename, pseudo_paths )
            if pseudo_filename is None: 
                raise ValueError('Pseudopotential %s not found in any of these paths:\n'%pseudo_filename+ppstring)

            if verbose: print('cp %s %s'%(pseudo_filename,destpath))
            shutil.copy(pseudo_filename,destpath)

    def set_kpoints(self,kpoints,shiftk=[0,0,0]):
        """Add some logic to set the kpoints mesh"""
        #sanity check
        if len(kpoints) != 3: raise ValueError('Wrong kpoints dimensions')
        self.kpoints = kpoints
        self.shiftk = shiftk

    def copy(self):
        """Return a copy of this instance"""
        import copy
        return copy.deepcopy(self)

    def read_atomicspecies(self):
        lines = iter(self.file_lines)
        #find ATOMIC_SPECIES keyword in file and read next line
        for line in lines:
            if "ATOMIC_SPECIES" in line:
                for i in range(int(self.system["ntyp"])):
                    atype, mass, psp = next(lines).split()
                    self.atypes[atype] = [mass,psp]

    def get_symmetry_spglib(self):
        """
        get the symmetry group of this system using spglib
        """
        import spglib

        lat, positions, atypes = self.get_atoms()
        lat = np.array(lat)

        at = np.unique(atypes)
        an = dict(list(zip(at,list(range(len(at))))))
        atypes = [an[a] for a in atypes]

        cell = (lat,positions,atypes)

        spacegroup = spglib.get_spacegroup(cell,symprec=1e-5)
        return spacegroup

    def get_masses(self):
        """ Get an array with the masses of all the atoms
        """
        masses = []
        for atom in self.atoms:
            atype = self.atypes[atom[0]]
            mass = float(atype[0])
            masses.append(mass)
        return masses

    def set_path(self,path):
        self.klist = path.get_klist()

    def get_atoms(self):
        """ Get the lattice parameters, postions of the atoms and chemical symbols
        """
        self.read_cell_parameters()
        cell = self.cell_parameters
        sym = [atom[0] for atom in self.atoms]
        pos = [atom[1] for atom in self.atoms]
        if self.atomic_pos_type == 'bohr':
            pos = car_red(pos,cell)
        return cell, pos, sym

    def set_atoms_string(self,string):
        """
        set the atomic postions using string of the form
        Si 0.0 0.0 0.0
        Si 0.5 0.5 0.5
        """
        atoms_str = [line.strip().split() for line in string.strip().split('\n')]
        self.atoms = []
        for atype,x,y,z in atoms_str:
            self.atoms.append([atype,list(map(float,[x,y,z]))])

    def set_atoms_ase(self,atoms):
        """ set the atomic postions using a Atoms datastructure from ase
        """
        # we will write down the cell parameters explicitly
        self.system['ibrav'] = 0
        if 'celldm(1)' in self.system: del self.system['celldm(1)']
        self.cell_parameters = atoms.get_cell()
        self.atoms = list(zip(atoms.get_chemical_symbols(),atoms.get_scaled_positions()))
        self.system['nat'] = len(self.atoms)

    def displace(self,mode,displacement,masses=None):
        """ A routine to displace the atoms acoording to a phonon mode
        """
        if masses is None:
            masses = [1] * len(self.atoms)
            small_mass = 1
        else:
            small_mass = min(masses) #we scale all the displacements to the bigger mass
        for i in range(len(self.atoms)):
            self.atoms[i][1] = self.atoms[i][1] + mode[i].real*displacement*sqrt(small_mass)/sqrt(masses[i])

    def read_atoms(self):
        lines = iter(self.file_lines)
        #find READ_ATOMS keyword in file and read next lines
        for line in lines:
            if "ATOMIC_POSITIONS" in line:
                atomic_pos_type = line
                self.atomic_pos_type = re.findall('([A-Za-z]+)',line)[-1]
                for i in range(int(self.system["nat"])):
                    atype, x,y,z = next(lines).split()
                    self.atoms.append([atype,[float(i) for i in (x,y,z)]])
        self.atomic_pos_type = atomic_pos_type.replace('{','').replace('}','').strip().split()[1]

    def read_cell_parameters(self):
        ibrav = int(self.system['ibrav'])
        def rmchar(string,symbols): return ''.join([c for c in string if c not in symbols])

        if ibrav == 0:
            if 'celldm(1)' in list(self.system.keys()):
                a = float(self.system['celldm(1)'])
            else:
                a = 1
            lines = iter(self.file_lines)
            for line in lines:
                if "CELL_PARAMETERS" in line:
                    units = rmchar(line.strip(),'{}()').split()
                    self.cell_parameters = [[],[],[]]
                    if len(units) > 1:
                        self.cell_units = units[1]
                    else:
                        self.cell_units = 'bohr'
                    for i in range(3):
                        self.cell_parameters[i] = [ float(x)*a for x in next(lines).split() ]
            if self.cell_units == 'angstrom' or self.cell_units == 'bohr':
                if 'celldm(1)' in self.system: del self.system['celldm(1)']
            if 'celldm(1)' not in list(self.system.keys()):
                a = np.linalg.norm(self.cell_parameters[0])
        elif ibrav == 1:
            a = float(self.system['celldm(1)'])
            self.cell_parameters = [[  a,   0,   0],
                                    [  0,   a,   0],
                                    [  0,   0,   a]]
        elif ibrav == 2:
            a = float(self.system['celldm(1)'])
            self.cell_parameters = [[ -a/2,   0, a/2],
                                    [    0, a/2, a/2],
                                    [ -a/2, a/2,   0]]
        elif ibrav == 3:
            a = float(self.system['celldm(1)'])
            self.cell_parameters = [[ a/2,  a/2,  a/2],
                                    [-a/2,  a/2,  a/2],
                                    [-a/2, -a/2,  a/2]]
        elif ibrav == 4:
            a = float(self.system['celldm(1)'])
            c = float(self.system['celldm(3)'])
            self.cell_parameters = [[   a,          0,  0],
                                    [-a/2,sqrt(3)/2*a,  0],
                                    [   0,          0,c*a]]
        elif ibrav == 6:
            a = float(self.system['celldm(1)'])
            c = float(self.system['celldm(3)'])
            self.cell_parameters = [[  a,   0,   0],
                                    [  0,   a,   0],
                                    [  0,   0, c*a]]
        else:
            raise NotImplementedError('ibrav = %d not implemented'%ibrav)
        self.alat = a 
        
    def read_kpoints(self):
        lines = iter(self.file_lines)
        #find K_POINTS keyword in file and read next line
        for line in lines:
            if "K_POINTS" in line:
                #chack if the type is automatic
                if "automatic" in line:
                    self.ktype = "automatic"
                    vals = list(map(float, next(lines).split()))
                    self.kpoints, self.shiftk = vals[0:3], vals[3:6]
                #otherwise read a list
                else:
                    #read number of kpoints
                    nkpoints = int(next(lines).split()[0])
                    self.klist = []
                    self.ktype = ""
                    try:
                        lines_list = list(lines)
                        for n in range(nkpoints):
                            vals = lines_list[n].split()[:4]
                            self.klist.append( list(map(float,vals)) )
                    except IndexError:
                        print("wrong k-points list format")
                        exit()

    def slicefile(self, keyword):
        lines = re.findall('&%s(?:.?)+\n((?:.+\n)+?)(?:\s+)?[\/&]'%keyword,"".join(self.file_lines),re.MULTILINE)
        return lines

    def store(self,group,name):
        """
        Save the variables specified in each of the groups on the structure
        """
        for file_slice in self.slicefile(name):
            for keyword, value in re.findall('([a-zA-Z_0-9_\(\)]+)(?:\s+)?=(?:\s+)?([a-zA-Z/\'"0-9_.-]+)',file_slice):
                group[keyword.strip()]=value.strip()

    def stringify_group(self, keyword, group):
        if group != {}:
            string='&%s\n' % keyword
            for keyword in sorted(group): # Py2/3 discrepancy in keyword order
                string += "%20s = %s\n" % (keyword, group[keyword])
            string += "/&end\n"
            return string
        else:
            return ''

    def remove_key(self,group,key):
        """ if a certain key exists in the group, remove it
        """
        if key in list(group.items()):
            del group[key]

    def run(self,filename,procs=1,folder='.'):
        """ this function is used to run this job locally
        """
        os.system('mkdir -p %s'%folder)
        self.write("%s/%s"%(folder,filename))
        if procs == 1:
            os.system('cd %s; OMP_NUM_THREADS=1 %s -inp %s > %s.log' % (folder,self._pw,filename,filename))
        else:
            os.system('cd %s; OMP_NUM_THREADS=1 mpirun -np %d %s -inp %s > %s.log' % (folder,procs,self._pw,filename,filename))

    def write(self,filename):
        f = open(filename,'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        """
        Output the file in the form of a string
        """
        string = ''
        string += self.stringify_group("control",self.control) #print control
        string += self.stringify_group("system",self.system) #print system
        string += self.stringify_group("electrons",self.electrons) #print electrons
        string += self.stringify_group("ions",self.ions) #print ions
        string += self.stringify_group("cell",self.cell) #print ions

        #print atomic species
        string += "ATOMIC_SPECIES\n"
        for atype in self.atypes:
            string += " %3s %8s %20s\n" % (atype, self.atypes[atype][0], self.atypes[atype][1])
        #print atomic positions
        string += "ATOMIC_POSITIONS { %s }\n"%self.atomic_pos_type
        for atom in self.atoms:
            string += "%3s %14.10lf %14.10lf %14.10lf\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2])
        #print kpoints
        if self.ktype == "automatic":
            string += "K_POINTS { %s }\n" % self.ktype
            string += ("%3d"*6+"\n")%tuple(self.kpoints + self.shiftk)
        elif self.ktype == "crystal":
            string += "K_POINTS { %s }\n" % self.ktype
            string += "%d\n" % len(self.klist)
            for i in self.klist:
              string += ('%12.8lf '*4+'\n') % tuple(i)
        else:
            string += "K_POINTS { }\n"
            string += "%d\n" % len(self.klist)
            for i in self.klist:
                string += (("%12.8lf "*4)+"\n")%tuple(i)
        if self.system['ibrav'] == 0 or self.system['ibrav'] == '0':
            string += "CELL_PARAMETERS %s\n"%self.cell_units
            for i in range(3):
                string += ("%14.10lf "*3+"\n")%tuple(self.cell_parameters[i])
        return string
