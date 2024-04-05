# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from qepy.auxiliary import *
from .lattice import *
from yambopy.plot.plotting import add_fig_kwargs 
from re import findall
from numpy import zeros

__all__ = ['PwXML']

HatoeV = 27.2107

class PwXML():
    # This class reads up to version 6.7
    """ Class to read data from a Quantum espresso XML file
    """
    _eig_xml   = 'eigenval.xml'
    _eig1_xml  = 'eigenval1.xml'
    _eig2_xml  = 'eigenval2.xml'

    
    def __init__(self,prefix,path='.',verbose=0):
        """ Initialize the structure with the path where the datafile.xml is
        """
        self.prefix = prefix
        self.path   = path

        datafiles = {'data-file.xml':        self.read_datafile,
                     'data-file-schema.xml': self.read_datafile_schema}

        done_reading = False

        #check if the name is data-file.xml or data-file-schema.xml or whatever....
        for filename,read in list(datafiles.items()):
            path_filename = "%s/%s.save/%s"%(path, prefix, filename)
            if os.path.isfile(path_filename):
                if verbose: print("reading %s"%filename)
                done_reading = read(path_filename)
                break

        #trap errors
        if not done_reading:
            possible_files = " or ".join(list(datafiles.keys()))
            raise ValueError('Failed to read %s in %s/%s.save'%(possible_files,path,prefix))

    def read_datafile(self,filename):
        """
        Read some data from the xml file in the old format of quantum espresso
        """
        self.datafile_xml = ET.parse( filename ).getroot()


        #get magnetization state

        self.lsda = False
        if 'T' in self.datafile_xml.findall("SPIN/LSDA")[0].text:
            self.lsda = True

        #get acell
        self.celldm = [ float(x) for x in self.datafile_xml.findall("CELL/CELL_DIMENSIONS")[0].text.strip().split('\n') ]

        #get cell
        self.cell = []
        for i in range(1,4):
            cell_lat = self.datafile_xml.findall("CELL/DIRECT_LATTICE_VECTORS/a%d"%i)[0].text
            self.cell.append([float(x) for x in cell_lat.strip().split()])

        #get reciprocal cell
        self.rcell = []
        for i in range(1,4):
            rcell_lat = self.datafile_xml.findall("CELL/RECIPROCAL_LATTICE_VECTORS/b%d"%i)[0].text
            self.rcell.append([float(x) for x in rcell_lat.strip().split()])

        #get atoms
        self.natoms = int(self.datafile_xml.findall("IONS/NUMBER_OF_ATOMS")[0].text)
        self.atoms = []
        for i in range(1,self.natoms+1):
            atom_xml = self.datafile_xml.findall("IONS/ATOM.%d"%i)[0]
            #read postions
            pos_string = atom_xml.get('tau').strip().split()
            pos = [float(x) for x in pos_string]
            #read type  
            atype = atom_xml.get('SPECIES').strip()
            #store
            self.atoms.append([atype,pos])

        #get atomic species
        self.natypes = int(self.datafile_xml.findall("IONS/NUMBER_OF_SPECIES")[0].text) 
        self.atypes = {}
        for i in range(1,self.natypes+1):
            atype_xml = self.datafile_xml.findall("IONS/SPECIE.%d"%i)[0]
            #read string
            atype_string = atype_xml.findall('ATOM_TYPE')[0].text.strip()
            #read mass
            atype_mass =  float(atype_xml.findall('MASS')[0].text.strip())
            #read pseudo
            atype_pseudo =  atype_xml.findall('PSEUDO')[0].text.strip()
            self.atypes[atype_string]=[atype_mass,atype_pseudo]

        #get nkpoints

        #get nkpoints
        self.nkpoints = int(self.datafile_xml.findall("BRILLOUIN_ZONE/NUMBER_OF_K-POINTS")[0].text.strip())
        # Read the number of BANDS
        self.nbands   = int(self.datafile_xml.find("BAND_STRUCTURE_INFO/NUMBER_OF_BANDS").text)

        #get k-points
        self.kpoints = [] 
        for i in range(self.nkpoints):
          k_aux = self.datafile_xml.findall('BRILLOUIN_ZONE/K-POINT.%d'%(i+1))[0].get('XYZ')
          self.kpoints.append([float(x) for x in k_aux.strip().split()])
        
        #get fermi
        self.fermi = float(self.datafile_xml.find("BAND_STRUCTURE_INFO/FERMI_ENERGY").text)*HatoeV
 
        #get eigenvalues

        if not self.lsda:

           eigen = []
           for ik in range(self.nkpoints):
               for EIGENVALUES in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig_xml) ).getroot().findall("EIGENVALUES"):
                   eigen.append(list(map(float,EIGENVALUES.text.split())*HatoeV))
           self.eigen  = eigen - self.fermi
           self.eigen1 = eigen - self.fermi

        #get eigenvalues of spin up & down

        if self.lsda:
           eigen1, eigen2 = [], []
           for ik in range(self.nkpoints):
               for EIGENVALUES1 in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig1_xml) ).getroot().findall("EIGENVALUES"):
                    eigen1.append(list(map(float, EIGENVALUES1.text.split())*HatoeV))
               for EIGENVALUES2 in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig2_xml) ).getroot().findall("EIGENVALUES"):
                    eigen2.append(list(map(float, EIGENVALUES2.text.split())*HatoeV))

           self.eigen   = eigen1 - self.fermi
           self.eigen1  = eigen1 - self.fermi
           self.eigen2  = eigen2 - self.fermi

        #get occupations of spin up & down
        if self.lsda:
           occ1, occ2 = [], []
           for ik in range(self.nkpoints):
               for OCCUPATIONS1 in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig1_xml) ).getroot().findall("OCCUPATIONS"):
                    occ1.append(list(map(float, OCCUPATIONS1.text.split())))
               for OCCUPATIONS2 in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig2_xml) ).getroot().findall("OCCUPATIONS"):
                    occ2.append(list(map(float, OCCUPATIONS2.text.split())))
        
           self.occupation1 = occ1
           self.occupation2 = occ2

        #get Bravais Lattice
        self.bravais_lattice = str(self.datafile_xml.find("CELL/BRAVAIS_LATTICE").text)
        if all(s in self.bravais_lattice for s in ["cubic","P"]):               self.ibrav = 1
        if all(s in self.bravais_lattice for s in ["cubic","F"]):               self.ibrav = 2
        if all(s in self.bravais_lattice for s in ["cubic","I"]):               self.ibrav = 3
        if all(s in self.bravais_lattice for s in ["Hexagonal","Trigonal"]):    self.ibrav = 4
        if all(s in self.bravais_lattice for s in ["Trigonal","R"]):            self.ibrav = 5
        if all(s in self.bravais_lattice for s in ["Tetragonal","P"]):          self.ibrav = 6
        if all(s in self.bravais_lattice for s in ["Tetragonal","I"]):          self.ibrav = 7
        if all(s in self.bravais_lattice for s in ["Orthorhombic","P"]):        self.ibrav = 8
        if all(s in self.bravais_lattice for s in ["Orthorhombic","base-centered"]): self.ibrav = 9
        if all(s in self.bravais_lattice for s in ["Orthorhombic","face-centered"]): self.ibrav = 10
        if all(s in self.bravais_lattice for s in ["Orthorhombic","body-centered"]): self.ibrav = 11
        if all(s in self.bravais_lattice for s in ["Monoclinic","P"]):          self.ibrav = 12
        if all(s in self.bravais_lattice for s in ["Monoclinic","base-centered"]):   self.ibrav = 13
        if all(s in self.bravais_lattice for s in ["Triclinic"]):               self.ibrav = 14


        return True

    def read_datafile_schema(self,filename):
        """
        Read the data from the xml file in the new format of quantum espresso
        """
        self.datafile_xml = ET.parse( filename ).getroot()

        # occupation type

        self.occ_type = self.datafile_xml.findall("input/bands/occupations")[0].text

        #get magnetization state
        # TO BE DONE!!!
        self.lsda = False
        if 'true' in self.datafile_xml.findall("input/spin/lsda")[0].text:
            self.lsda = True
            #raise ValueError('Spin states not yet implemented for data-file-schema.xml') 

        #get cell
        self.cell = []
        for i in range(1,4):
            cell_lat = self.datafile_xml.findall("output/atomic_structure/cell/a%d"%i)[0].text
            self.cell.append([float(x) for x in cell_lat.strip().split()])

        #calculate acell
        self.acell = [ np.linalg.norm(a) for a in self.cell ]

        #get reciprocal cell
        self.rcell = []
        for i in range(1,4):
            rcell_lat = self.datafile_xml.findall("output/basis_set/reciprocal_lattice/b%d"%i)[0].text
            self.rcell.append([float(x) for x in rcell_lat.strip().split()])

        #get atoms
        self.natoms = int(self.datafile_xml.findall("output/atomic_structure")[0].get('nat'))
        self.atoms = []
        atoms = self.datafile_xml.findall("output/atomic_structure/atomic_positions/atom")
        for i in range(self.natoms):
            atype = atoms[i].get('name')
            pos = [float(x) for x in atoms[i].text.strip().split()]
            self.atoms.append([atype,pos])

        #get atomic species
        self.natypes = int(self.datafile_xml.findall("output/atomic_species")[0].get('ntyp')) 
        atypes = self.datafile_xml.findall("output/atomic_species/species")
        self.atypes = {}
        for i in range(self.natypes):
            #read string
            atype_string = atypes[i].get('name').strip()
            #read mass
            atype_mass = atypes[i].findall('mass')[0].text.strip()
            #read pseudo
            atype_pseudo = atypes[i].findall('pseudo_file')[0].text.strip()
            self.atypes[atype_string]=[atype_mass,atype_pseudo]

        #get nkpoints
        self.nkpoints = int(self.datafile_xml.findall("output/band_structure/nks")[0].text.strip())
        # Read the number of BANDS
        if self.lsda:
           self.nbands_up = int(self.datafile_xml.findall("output/band_structure/nbnd_up")[0].text.strip())
           self.nbands_dw = int(self.datafile_xml.findall("output/band_structure/nbnd_dw")[0].text.strip())
           self.nbands = self.nbands_up + self.nbands_dw
        else:
           self.nbands = int(self.datafile_xml.findall("output/band_structure/nbnd")[0].text.strip())

        #get ks states
        kstates = self.datafile_xml.findall('output/band_structure/ks_energies')

        #get k-points
        self.kpoints = [] 
        for i in range(self.nkpoints):
            kpoint = [float(x) for x in kstates[i].findall('k_point')[0].text.strip().split()]
            self.kpoints.append( kpoint )

        #get fermi (it depends on the occupations)
        if self.occ_type == 'fixed':
           self.fermi = float(self.datafile_xml.find("output/band_structure/highestOccupiedLevel").text)*HatoeV
        else:
           self.fermi = float(self.datafile_xml.find("output/band_structure/fermi_energy").text)*HatoeV

        #get eigenvalues
        self.eigen1 = []
        for k in range(self.nkpoints):
            eigen = [float(x) for x in kstates[k].findall('eigenvalues')[0].text.strip().split()]
            self.eigen1.append( eigen )
        self.eigen1 = np.array(self.eigen1)*HatoeV - self.fermi
 
        #get Bravais lattice
        self.ibrav = self.datafile_xml.findall("output/atomic_structure")[0].get('bravais_index')

        return True

    def get_scaled_atoms(self):
        """Build and return atoms dictionary"""
        atoms = []
        for atype,pos in self.atoms:
            red_pos = car_red([pos],self.cell)[0].tolist()
            atoms.append([atype,red_pos])
        return atoms
        
    def get_atypes_dict(self):
        return self.atypes

    def get_lattice_dict(self):
        return dict(ibrav=self.ibrav,cell_parameters=self.cell)

    def get_structure_dict(self):
        """
        Get a structure dict that can be used to create a new pw input file
        """
        return dict(atoms=self.get_scaled_atoms(),
             atypes=self.get_atypes_dict(),
             lattice=self.get_lattice_dict())

    def get_cartesian_positions(self):
        """ get the atomic positions in cartesian coordinates """
        positions = []
        for atype,pos in self.atoms:
            positions.append(pos)
        return np.array(positions)

    def get_scaled_positions(self):
        """ get the atomic positions in reduced coordinates """
        cartesian_positions = self.get_cartesian_positions()
        return car_red(cartesian_positions,self.cell)

    def __str__(self):
        lines = []; app = lines.append
        app("cell:")
        for c in self.cell:
            app(("%12.8lf "*3)%tuple(c))
        app("atoms:")
        for atype,pos in self.atoms:
            app("%3s"%atype + ("%12.8lf "*3)%tuple(pos))
        app("nkpoints: %d"%self.nkpoints)
        app("nbands:   %d"%self.nbands)
        return "\n".join(lines)

    def plot_eigen_ax(self,ax,path_kpoints=[],xlim=(),ylim=(),color='r',**kwargs):
        #
        # Careful with variable path. I am substituting vy path_kpoints
        # To be done in all the code (and in the tutorials)
        #
        # argurments:
        # ls: linestyle
        if path_kpoints:
            if isinstance(path_kpoints,Path):
                path_kpoints = path_kpoints.get_indexes()
                path_ticks, path_labels = list(zip(*path_kpoints))
            ax.set_xticks( path_ticks )
            ax.set_xticklabels( path_labels )
        ax.set_ylabel('E (eV)')

        ls = kwargs.pop('ls','solid')
        lw = kwargs.pop('lw',1)
        y_offset = kwargs.pop('y_offset',0.0)
        #get kpoint_dists 
        kpoints_dists = calculate_distances(self.kpoints)
        ticks, labels = list(zip(*path_kpoints))
        ax.set_xticks([kpoints_dists[t] for t in ticks])
        ax.set_xticklabels(labels)
        ax.set_xlim(kpoints_dists[0],kpoints_dists[-1])

        #plot vertical lines
        for t in ticks:
            ax.axvline(kpoints_dists[t],c='k',lw=2)
        ax.axhline(0,c='k')

        #plot bands
       
        if self.lsda:
           eigen1 = np.array(self.eigen1)

           for ib in range(self.nbands_up):
               ax.plot(kpoints_dists,eigen1[:,ib]                + y_offset, '%s-'%color, lw=lw, zorder=1,label='spin-up') # spin-up 
               ax.plot(kpoints_dists,eigen1[:,ib+self.nbands_up] + y_offset, 'b-', lw=lw, zorder=1,label='spin-down') # spin-down

           import matplotlib.pyplot as plt
           handles, labels = plt.gca().get_legend_handles_labels()
           by_label = dict(zip(labels, handles))
           plt.legend(by_label.values(), by_label.keys())

        # Case: Non spin polarization
        else:
           eigen1 = np.array(self.eigen1)

           for ib in range(self.nbands):
               ax.plot(kpoints_dists,eigen1[:,ib] + y_offset, color=color,linestyle=ls , lw=lw, zorder =1)

        #plot options
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)

     
    #def plot_eigen_spin_ax(self,ax,path_kpoints=[],xlim=(),ylim=(),spin_proj=None):
    def plot_eigen_spin_ax(self,ax,path_kpoints=[],xlim=(),ylim=(),spin_proj=False,spin_folder='.'):
        #
        # Careful with variable path. I am substituting vy path_kpoints
        # To be done in all the code (and in the tutorials)
        # This is a test function for spin-polarized bands
        #
        #
        
        import matplotlib.pyplot as plt
        self.spin_proj = np.array(spin_proj) if spin_proj is not None else None 

        if spin_proj == True:
           self.spin_projection(spin_dir=3,folder=spin_folder,prefix='bands')
           print(self.spin_3)
        exit()

        if path_kpoints:
            if isinstance(path_kpoints,Path):
                path_kpoints = path_kpoints.get_indexes()
                path_ticks, path_labels = list(zip(*path_kpoints))
            ax.set_xticks( path_ticks )
            ax.set_xticklabels( path_labels )
        ax.set_ylabel('E (eV)')

        # I choose a colormap for spin
        color_map  = plt.get_cmap('seismic')

        #get kpoint_dists 
        kpoints_dists = calculate_distances(self.kpoints)
        ticks, labels = list(zip(*path_kpoints))
        ax.set_xticks([kpoints_dists[t] for t in ticks])
        ax.set_xticklabels(labels)
        ax.set_xlim(kpoints_dists[0],kpoints_dists[-1])

        # NOT WORKING, CHECK IT!
        #plot vertical lines
        #for t in ticks:
        #    ax.axvline(kpoints_dists[t],c='k',lw=2)
        #ax.axhline(0,c='k')

        #plot bands
        eigen1 = np.array(self.eigen1)
        for ib in range(self.nbands):
            x = kpoints_dists
            y = eigen1[:,ib] - self.fermi
            color_spin = self.spin_proj[:,ib] + 0.5 # I renormalize 0 => down; 1 => up
            ax.scatter(x,y,s=100,c=color_spin,cmap=color_map,vmin=0.0,vmax=1.0,edgecolors='none')
       
        #plot spin-polarized bands: TO BE DONE
        #if self.lsda:

        #  eigen2 = np.array(self.eigen2)
        #   for ib in range(self.nbands):
        #       ax.plot(kpoints_dists,eigen2[:,ib]*HatoeV - self.fermi*HatoeV, 'b-', lw=2)

        #plot options
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)


    '''
    Workaround to include occupations in the plot. AMS
    '''

    def plot_eigen_occ_ax(self,ax,path_kpoints=[],xlim=(),ylim=(),color='r'):

        if path_kpoints:
            if isinstance(path_kpoints,Path):
                path_kpoints = path_kpoints.get_indexes()
            ax.set_xticks( *list(zip(*path_kpoints)) )
        ax.set_ylabel('E (eV)')

        #get kpoint_dists 
        kpoints_dists = calculate_distances(self.kpoints)
        ticks, labels = list(zip(*path_kpoints))
        ax.set_xticks([kpoints_dists[t] for t in ticks])
        ax.set_xticklabels(labels)
        ax.set_xlim(kpoints_dists[0],kpoints_dists[-1])

        #plot vertical lines
        for t in ticks:
            ax.axvline(kpoints_dists[t],c='k',lw=2)
        ax.axhline(0,c='k')
        import matplotlib.pyplot as plt

        #plot bands
        eigen1 = np.array(self.eigen1)
        occ1   = np.array(self.occupation1)
        for ib in range(self.nbands):
            plt.scatter(kpoints_dists,eigen1[:,ib] - self.fermi, s=10*occ1[:,ib],c=color)
       
        #plot spin-polarized bands
        if self.lsda:

           eigen2 = np.array(self.eigen2)
           occ2   = np.array(self.occupation1)
           for ib in range(self.nbands):
               plt.scatter(kpoints_dists,eigen2[:,ib] - self.fermi, s=10*occ2[:,ib],c='b')


        #plot options
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)

    @add_fig_kwargs
    def plot_eigen(self,path_kpoints=[],xlim=(),ylim=()):
        """ plot the eigenvalues using matplotlib
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_eigen_ax(ax,path_kpoints=path_kpoints)
        return fig

    def write_eigen(self,fmt='gnuplot'):
        """ write eigenvalues to a text file
        """
        if fmt=='gnuplot':
            f = open('%s.dat'%self.prefix,'w')
            for ib in range(self.nbands):
                for ik in range(self.nkpoints):
                    f.write("%.1lf %.4lf \n " % (ik,self.eigen1[ik][ib]) )
                f.write("\n")
            f.close()
        else:
            print('fmt %s not implemented'%fmt)

    def spin_projection(self,spin_dir=3,folder='.',prefix='bands'):
        """
        This function reads the spin projection given by bands.x in txt file
        lsigma(i) = .true.
        By default I set the spin direction z ==3

        """
        if spin_dir ==3:  
        #data_eigen  = open('%s/%s.out'  % (folder,prefix),'r').readlines()
           data_spin_3 = open('%s/%s.out.3'% (folder,prefix),'r').readlines()
           
           # check consistency file from bands.x and xml file
           nband = int(findall(r"[-+]?\d*\.\d+|\d+", data_spin_3[0].strip().split()[2]  )[0])
           nk    = int(data_spin_3[0].strip().split()[-2])
           nline = int(nband/10)
           if nband < 10: print("Error, uses only nband => 10 and multiple of 10")
           if self.nbands != nband or self.nkpoints != nk: print("Warning: Dimensions are different!")

           self.spin_3 = zeros([self.nkpoints,self.nbands])

        for ik in range(self.nkpoints):
            for ib in range(nline):
                ib1, ib2, ib3 = int(ib*10), int((ib+1)*10), int(ik*(nband/10+1)+2+ib)
                self.spin_3[ik,ib1:ib2] = list( map(float,data_spin_3[ib3].split()))

    def read_symmetries(self):
        """
        Read symmetry operations from data-file-schema.xml

        Works for ibrav>0

        NB: data-file-schema.xml has nrot and nsym with nsym<nrot.
            They are stored in order with nsym first.

        NB2: Most likely not working with symmorphic symmetries
        """
        #nsym
        nrot = int(self.datafile_xml.findall("output/symmetries/nrot")[0].text.strip())        
        self.nsym = int(self.datafile_xml.findall("output/symmetries/nsym")[0].text.strip())        
        no_t_rev = ( self.datafile_xml.findall("input/symmetry_flags/no_t_rev")[0].text.strip() == "true" )
        self.nsym_with_trev = 2*self.nsym 

        # data
        if not no_t_rev: sym_red = np.zeros((self.nsym_with_trev,3,3))
        else:            sym_red = np.zeros((self.nsym,3,3))
        symmetries = self.datafile_xml.findall("output/symmetries/symmetry/rotation") # All rotations
        symmetries = symmetries[:self.nsym] # Retain only point group symmetries with no trev
        #sym_info   = self.datafile_xml.findall("output/symmetries/symmetry/info")
        #NB: sym_info[:].attrib['class'] contains irrep names if found

        #read (non-symmorphic) symmetries
        for i in range(self.nsym):
            sym = np.array( [float(x) for x in symmetries[i].text.strip().split()] ).reshape(3,3)
            sym_red[i] = sym.T # symmetries are saved as the transposed in the .xml 
        self.sym_red = sym_red.astype(int)

        #convert to c.c.
        self.sym_car = self.sym_red_car()

        #check for time reversal and apply it to sym_car
        if not no_t_rev: self.apply_t_rev()

    def apply_t_rev(self):
        """
        Add T*S rotation matrices
        """
        for n in range(self.nsym,self.nsym_with_trev): self.sym_car[n]=-1*self.sym_car[n-self.nsym]

    def sym_red_car(self):
        """
        Transform symmetry ops. in Cartesian coordinates
        """
        lat = np.array(self.cell)
        sym_car = np.zeros([len(self.sym_red),3,3],dtype=float)
        for n,s in enumerate(self.sym_red):
            sym_car[n] = np.dot( np.linalg.inv(lat), np.dot(s, lat ) ).T
        return sym_car

    def expand_kpoints_xml(self,atol=1e-6,expand_eigen=True,verbose=0):
        """
        Take a list of kpoints and symmetry operations and return the full brillouin zone
        with the corresponding index in the irreducible brillouin zone

        Expand also eigenvalues by default
        """
        alat      = np.array(self.acell)
        kpts_ibz  = np.array(self.kpoints)
        eigen_ibz = np.array(self.eigen1)
        rlat      = np.array(self.rcell)

        self.read_symmetries()

        #check if the kpoints were already exapnded
        kpoints_indices  = []
        kpoints_full     = []
        symmetry_indices = []

        #kpoints in the full brillouin zone organized per index
        kpoints_full_i = {}

        #expand using symmetries
        for nk,k in enumerate(kpts_ibz):
            #if the index in not in the dictionary add a list
            if nk not in kpoints_full_i:
                kpoints_full_i[nk] = []

            for ns,sym in enumerate(self.sym_car):

                new_k = np.dot(sym,k)

                #check if the point is inside the bounds
                k_red = car_red([new_k],rlat)[0]
                k_bz = (k_red+atol)%1

                #if the vector is not in the list of this index add it
                if not vec_in_list(k_bz,kpoints_full_i[nk]):
                    kpoints_full_i[nk].append(k_bz)
                    kpoints_full.append(new_k)
                    kpoints_indices.append(nk)
                    symmetry_indices.append(ns)
                    continue

        #calculate the weights of each of the kpoints in the irreducible brillouin zone
        nkpoints_full = len(kpoints_full)
        weights = np.zeros([nkpoints_full])
        for nk in kpoints_full_i:
            weights[nk] = float(len(kpoints_full_i[nk]))/nkpoints_full

        if verbose: print("%d kpoints expanded to %d"%(self.nkpoints,len(kpoints_full)))

        #set the variables
        self.weights_ibz      = np.array(weights)
        self.kpoints_indices  = np.array(kpoints_indices)
        self.symmetry_indices = np.array(symmetry_indices)
        self.nkbz             = nkpoints_full
        #cartesian coordinates of QE
        self.kpoints_bz       = np.array(kpoints_full)

        if expand_eigen:

            self.eigen_bz = np.zeros((self.nkbz,self.nbands))
            for ik in range(self.nkbz): self.eigen_bz[ik,:] = eigen_ibz[self.kpoints_indices[ik],:]
            if verbose: print("Eigenvalues also expanded.")
