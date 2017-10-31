#
# Copyright (C) 2017 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
import re
import numpy as np
from collections import defaultdict
from netCDF4 import Dataset
from .yambofile import YamboFile
from yambopy.jsonencoder import JsonDumper
from yambopy import YamboIn
from yambopy.units import ha2ev

class YamboOut():
    """ 
    Class to read yambo output files and pack them in a .json file

    **Arguments:**

    ``folder``:      The relative path of the folder where yambo dumped its input files

    ``save_folder``: The path were the SAVE folder is localized 
    """

    _tags = ['refl', 'eel', 'eps', 'qp', 'sf',
             'carriers', 'polarization', 'external']
    _netcdf = ['ndb.QP', 'ndb.HF_and_locXC']
    _tagsexp = r'#\n#\s+((?:(?:[`0-9a-zA-Z\-\/\|\\(\)\_[\]]+)\s+)+)#\n\s'

    def __init__(self, folder, save_folder='.'):
        """
        For now we initialize the class in the same way as before (from_files)
        In the future we should be able to initialize this class from
        a file in the disk or from a dictionary
        """
        self.files = defaultdict(dict)
        self.from_folder(folder,save_folder)

    def from_folder(self,folder,save_folder):
        self.folder = folder

        # check if the save folder is in save_folder if not try folder
        if not os.path.isdir(save_folder + '/SAVE'):
            self.save_folder = folder
        else:
            self.save_folder = save_folder

        # get all the files in the output dir
        if os.path.isdir(folder):
            outdir = os.listdir(folder)
        else:
            raise ValueError("Invalid folder: %s" % folder)

        # get the log dir
        logdir_path = os.path.join(folder, "LOG")
        if os.path.isdir(logdir_path):
            logdir = os.listdir(logdir_path)
        else:
            logdir = outdir

        # get output filenames
        self.netcdf = ["%s" % f for f in outdir if f in self._netcdf]
        self.output = ["%s" % f for f in outdir if f.startswith('o-') and YamboFile.has_tag(f,self._tags)]
        self.run = ["%s" % f for f in outdir if f.startswith('r-')]
        self.logs = ["%s" % f for f in logdir if f.startswith('l-')]

        # get data from output file
        self.get_runtime()
        self.get_outputfile()
        self.get_netcdffile()
        self.get_inputfile()
        self.get_cell()

    @staticmethod
    def has_output(folder):
        """Check if the folder has output files"""
        return [filename for filename in os.listdir(folder) if YamboFile.is_output(filename)] 

    def get_cell(self):
        """ 
        Get information about the unit cell (lattice vectors, atom types, positions,
        kpoints and symmetry operations) from the SAVE folder.
        """
        path = self.save_folder + '/SAVE/ns.db1'
        if os.path.isfile(path):
            self.nc_db = Dataset(path)
            ncvars = self.nc_db.variables

            self.lat = ncvars['LATTICE_VECTORS'][:].T
            self.alat = ncvars['LATTICE_PARAMETER'][:].T
            self.sym_car = ncvars['SYMMETRY'][:]
            self.kpts_iku = ncvars['K-POINTS'][:].T

            # read atoms
            atomic_pos = ncvars['ATOM_POS'][:]
            atomic_number = ncvars['atomic_numbers'][:].T.astype(int)
            self.nspecies = ncvars['number_of_atom_species'][0].astype(int)

            atom_positions = []
            atom_number = []
            for specie in range(self.nspecies):
                atoms_specie = atomic_pos[specie]
                for atom in atoms_specie:
                    atom_positions.append(atom)
                    atom_number.append(atomic_number[specie])

            self.atomic_positions = atom_positions
            self.atomic_number = atom_number
        else:
            raise ValueError('Could not find ns.db1 in %s from YamboOut' % path)

    def get_outputfile(self):
        """ 
        Get the data from the o-* files
        """
        # for all the o-* files
        for filename in self.output:
            
            #yambofile read
            # TODO: use Yambofile class to read the o-* files
            yf = YamboFile(filename,self.folder)

            #classic read
            with open(os.path.join(self.folder, filename),'r') as f:
                string = f.read()

                # get tags
                find_tags = re.findall(self._tagsexp, string)
                tags = [tag.strip() for tag in find_tags[0].split()]
                f.seek(0)

                # get data
                data = np.loadtxt(f, unpack=True)

                #store data
                self.files[filename].update(dict(zip(tags,data)))
                self.files[filename]["type"] = yf.type

    def get_netcdffile(self):
        """
        Get the netcdf files
            The supported netcdf files so far are:
                ndb.QP
                ndb.HF_and_locXC
        """
        for filename in self.netcdf:
            yf = YamboFile(filename, self.folder)
            #convert units
            if yf.type == 'netcdf_gw':
                yf.data['E-Eo'] *= ha2ev
                yf.data['Eo'] *= ha2ev
                yf.data['E'] *= ha2ev
            #save data
            self.files[filename] = yf.data
            self.files[filename]["type"] = yf.type

    def get_inputfile(self):
        """
        Get the input file from the o-* file
        """

        for filename in self.output:
            inputfile = []
            # read this inputfile
            with open(os.path.join(self.folder, filename),'r') as f:
                for line in f:
                    if 'Input file :' in line:
                        for line in f:
                            # Note this: to read the input file we just ignore the 
                            # first 4 characters of the section after the tag 'Input file:'
                            inputfile.append(line[4:])

            # use YamboIn to read the input file to a list
            yi = YamboIn(filename=None)
            self.files[filename]["input"] = yi.read_string(''.join(inputfile))

    def get_runtime(self):
        """
        Get the runtime from the r-* file
        """

        for filename in self.run:
            timing = {}

            with open(os.path.join(self.folder,filename),'r') as f:

                category = "UNKNOWN"
                for line in f:
                    if 'Timing' in line:
                        timing[category] = line.split()[-1].split('/')
                    if re.search('(\[[0-9.]+\].[A-Z])', line):
                        category = line.strip()
                
            self.files[filename]["runtime"] = timing
            self.files[filename]["type"] = "report"

    def get_data(self, tags):
        """
        Search for a tag in the output files and obtain the data
        """
        data = {}
        for key in self.data.keys():
            if all(tag in key for tag in tags):
                data[key] = dict(
                    list(zip(self.tags[key], np.array(self.data[key]).T)))
        return data

    def print_runtime(self):
        """
        Print the runtime in a string
        """
        timing = self.get_runtime()
        for t in list(timing.items()):
            print(t[0], '\n', t[1], '\n')

    def pack(self, filename=None):
        """
        Pack up all the data in the structure in a json file
        """
        # if no filename is specified we use the same name as the folder
        if not filename:
            filename = self.folder

        # create json dictionary
        jsondata = {"files": self.files,
                    "lattice": self.lat,
                    "alat": self.alat,
                    "kpts_iku": self.kpts_iku,
                    "sym_car": self.sym_car,
                    "atompos": self.atomic_positions,
                    "atomtype": self.atomic_number}

        filename = '%s.json' % filename
        JsonDumper(jsondata, filename)

    def __str__(self):
        lines = []
        lines.append("logs:")
        lines += self.logs
        lines.append("\nrun:")
        lines += self.run
        lines.append("\noutput (text):")
        lines += self.output
        lines.append("\noutput (netcdf):")
        lines += self.netcdf
        lines.append("\nlattice:")
        lines += [("%12.8lf " * 3) % tuple(vec) for vec in self.lat]
        lines.append("\natom positions:")
        for an, pos in zip(self.atomic_number, self.atomic_positions):
            lines.append( "%3d " % an + ("%12.8lf " * 3) % tuple(pos) )
        return "\n".join(lines)
