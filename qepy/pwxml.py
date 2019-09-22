# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
import xml.etree.ElementTree as ET
from qepy.auxiliary import *
from .lattice import *
from yambopy.plot.plotting import add_fig_kwargs 

__all__ = ['PwXML']

HatoeV = 27.2107

class PwXML():
    """ Class to read data from a Quantum espresso XML file
    """
    _eig_xml   = 'eigenval.xml'
    _eig1_xml  = 'eigenval1.xml'
    _eig2_xml  = 'eigenval2.xml'

    
    def __init__(self,prefix,path='.',verbose=0):
        """ Initlize the structure with the path where the datafile.xml is
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
        self.nkpoints = int(self.datafile_xml.findall("BRILLOUIN_ZONE/NUMBER_OF_K-POINTS")[0].text.strip())
        # Read the number of BANDS
        self.nbands   = int(self.datafile_xml.find("BAND_STRUCTURE_INFO/NUMBER_OF_BANDS").text)

        #get k-points
        self.kpoints = [] 
        for i in range(self.nkpoints):
          k_aux = self.datafile_xml.findall('BRILLOUIN_ZONE/K-POINT.%d'%(i+1))[0].get('XYZ')
          self.kpoints.append([float(x) for x in k_aux.strip().split()])
 
        #get eigenvalues

        if not self.lsda:

           eigen = []
           for ik in range(self.nkpoints):
               for EIGENVALUES in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig_xml) ).getroot().findall("EIGENVALUES"):
                   eigen.append(list(map(float, EIGENVALUES.text.split())))
           self.eigen  = eigen
           self.eigen1 = eigen

        #get eigenvalues of spin up & down

        if self.lsda:
           eigen1, eigen2 = [], []
           for ik in range(self.nkpoints):
               for EIGENVALUES in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig1_xml) ).getroot().findall("EIGENVALUES"):
                    eigen1.append(list(map(float, EIGENVALUES.text.split())))
               for EIGENVALUES in ET.parse( "%s/%s.save/K%05d/%s" % (self.path,self.prefix,(ik + 1),self._eig2_xml) ).getroot().findall("EIGENVALUES"):
                    eigen2.append(list(map(float, EIGENVALUES.text.split())))
           self.eigen1  = eigen1
           self.eigen2  = eigen2

        #get fermi
        self.fermi = float(self.datafile_xml.find("BAND_STRUCTURE_INFO/FERMI_ENERGY").text)

        #get Bravais Lattice
        # Just do for ibrav 4, should we do more general for all ibravs?
        self.bravais_lattice = str(self.datafile_xml.find("CELL/BRAVAIS_LATTICE").text)
        if all(s in self.bravais_lattice for s in ["cubic","P"]):            self.ibrav = 1
        if all(s in self.bravais_lattice for s in ["Hexagonal","Trigonal"]): self.ibrav = 4

        return True

    def read_datafile_schema(self,filename):
        """
        Read the data from the xml file in the new format of quantum espresso
        """
        self.datafile_xml = ET.parse( filename ).getroot()

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
        self.nbands = int(self.datafile_xml.findall("output/band_structure/nbnd")[0].text.strip())

        #get ks states
        kstates = self.datafile_xml.findall('output/band_structure/ks_energies')

        #get k-points
        self.kpoints = [] 
        for i in range(self.nkpoints):
            kpoint = [float(x) for x in kstates[i].findall('k_point')[0].text.strip().split()]
            self.kpoints.append( kpoint )

        #get eigenvalues
        self.eigen = []
        for k in range(self.nkpoints):
            eigen = [float(x) for x in kstates[k].findall('eigenvalues')[0].text.strip().split()]
            self.eigen.append( eigen )
        self.eigen = np.array(self.eigen)
 
        #get fermi
        self.fermi = float(self.datafile_xml.find("output/band_structure/highestOccupiedLevel").text)

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
        return dict(ibrav=0,cell_parameters=self.cell)

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

    def plot_eigen_ax(self,ax,path=[],xlim=(),ylim=()):
        print(path)
        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            ax.set_xticks( *list(zip(*path)) )
        ax.set_ylabel('E (eV)')

        #get kpoint_dists 
        kpoints_dists = calculate_distances(self.kpoints)
        ticks, labels = list(zip(*path))
        ax.set_xticks([kpoints_dists[t] for t in ticks])
        ax.set_xticklabels(labels)
        ax.set_xlim(kpoints_dists[0],kpoints_dists[-1])

        #plot vertical line
        #for point in path:
        #    x, label = point
        #    ax.axvline(x)
        #ax.axhline(0)

        #plot vertical lines
        for t in ticks:
            ax.axvline(kpoints_dists[t],c='k',lw=2)
        ax.axhline(0,c='k')

        #plot bands
        eigen = np.array(self.eigen1)
        for ib in range(self.nbands):
            ax.plot(kpoints_dists,eigen[:,ib]*HatoeV - self.fermi*HatoeV, 'r-', lw=2)
       
        #plot spin-polarized bands
        if self.lsda:

           eigen2 = np.array(self.eigen2)
           for ib in range(self.nbands):
               ax.plot(kpoints_dists,eigen2[:,ib]*HatoeV - self.fermi*HatoeV, 'b-', lw=2)


        #plot options
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)

    @add_fig_kwargs
    def plot_eigen(self,path=[],xlim=(),ylim=()):
        """ plot the eigenvalues using matplotlib
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_eigen_ax(ax,path=path)
        return fig

    def write_eigen(self,fmt='gnuplot'):
        """ write eigenvalues to a text file
        """
        if fmt=='gnuplot':
            f = open('%s.dat'%self.prefix,'w')
            for ib in range(self.nbands):
                for ik in range(self.nkpoints):
                    f.write("%.1lf %.4lf \n " % (ik,self.eigen[ik][ib]*HatoeV) )
                f.write("\n")
            f.close()
        else:
            print('fmt %s not implemented'%fmt)
