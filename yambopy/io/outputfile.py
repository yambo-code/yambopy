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
from netCDF4 import Dataset
from yamboparser import YamboFile
from yambopy.jsonencoder import JsonDumper
from yambopy import YamboIn

class YamboOut():
    """ 
    Class to read yambo output files and pack them in a .json file

    **Arguments:**

    ``folder``:      The relative path of the folder where yambo dumped its input files

    ``save_folder``: The path were the SAVE folder is localized 
    """

    _tags = ['refl','eel','eps','qp','sf','carriers','polarization','external']
    _netcdf = ['ndb.QP','ndb.HF_and_locXC']
    _tagsexp = r'#\n#\s+((?:(?:[`0-9a-zA-Z\-\/\|\\(\)\_[\]]+)\s+)+)#\n\s'

    def __init__(self,folder,save_folder='.'):

        self.folder = folder

        #check if the save folder is in save_folder if not try folder
        if not os.path.isdir(save_folder+'/SAVE'):
            self.save_folder = folder
        else:
            self.save_folder = save_folder

        #get all the files in the output dir
        if os.path.isdir(folder):
            outdir = os.listdir(folder)
        else:
            raise ValueError( "Invalid folder: %s"%folder )

        #get the log dir
        logdir_path = os.path.join(folder,"LOG")
        if os.path.isdir(logdir_path):
            logdir = os.listdir(logdir_path)
        else:
            logdir = outdir

        def has_tag(filename):
            """check if the filename has a tag in its name"""
            return any([tag in filename for tag in self._tags])

        #get output filenames 
        self.netcdf = ["%s"%f for f in outdir if f in self._netcdf ]
        self.output = ["%s"%f for f in outdir if f.startswith('o-') and has_tag(f)]
        self.run    = ["%s"%f for f in outdir if f.startswith('r-')]
        self.logs   = ["%s"%f for f in logdir if f.startswith('l-')]

        #get data from output file
        self.get_runtime()
        self.get_outputfile()
        self.get_netcdffile()
        self.get_inputfile()
        self.get_cell()
    
    def get_cell(self):
        """ 
        Get information about the unit cell (lattice vectors, atom types, positions,
        kpoints and symmetry operations) from the SAVE folder.
        """
        path = self.save_folder+'/SAVE/ns.db1'
        if os.path.isfile(path):
            self.nc_db         = Dataset(path)
            self.lat           = self.nc_db.variables['LATTICE_VECTORS'][:].T
            self.alat          = self.nc_db.variables['LATTICE_PARAMETER'][:].T
            self.sym_car       = self.nc_db.variables['SYMMETRY'][:]
            self.kpts_iku      = self.nc_db.variables['K-POINTS'][:].T
            
            #read atoms
            atomic_pos         = self.nc_db.variables['ATOM_POS'][:]
            atomic_number      = self.nc_db.variables['atomic_numbers'][:].T
            self.nspecies      = self.nc_db.variables['number_of_atom_species'][0].astype(int)

            atom_positions = []
            atom_number    = []
            for specie in range(self.nspecies):
                atoms_specie = atomic_pos[specie]
                for atom in atoms_specie:
                    atom_positions.append(atom)
                    atom_number.append(atomic_number[specie])

            self.atomic_positions = atom_positions
            self.atomic_number = atom_number
        else:
            raise ValueError('Could not find ns.db1 in %s from YamboOut'%path)

    def get_outputfile(self):
        """ 
        Get the data from the o-* files
        """

        self.data = {}
        self.tags = {}

        #for all the o-* files
        for filename in self.output:

            with open("%s/%s"%(self.folder,filename)) as f:
                string = f.read()

                #get tags
                tags = [tag.strip() for tag in re.findall(self._tagsexp,string)[0].split()]
                f.seek(0)
                self.tags[filename] = tags

                #get data
                self.data[filename] = np.loadtxt(f)

    def get_netcdffile(self):
        """
        Get the netcdf files
            The supported netcdf files so far are:
                ndb.QP
                ndb.HF_and_locXC
        """
        for filename in self.netcdf:
            yf = YamboFile(filename,self.folder)
            
            #from dictionary to tuple list
            tags, data = zip(*yf.data.items())

            self.data[filename] = data
            self.tags[filename] = tags

    def get_inputfile(self):
        """
        Get the input file from the o-* file
        """
        files = [open("%s/%s"%(self.folder,f),'r') for f in self.output]

        self.inputfile = {}
    
        for filename,f in zip(self.output,files):

            #read this inputfile
            inputfile = []
            for line in f:
                if 'Input file :' in line:
                    for line in f:
                        # Note this: to read the input file we just ignore the first 4 characters
                        # of the section after the tag 'Input file:'
                        inputfile.append( line[4:] )
            
            #use YamboIn to read the input file to a list
            yi = YamboIn(filename=None)
            self.inputfile[filename] = yi.read_string( ''.join(inputfile) )

        #close all the files
        for f in files: f.close()

    def get_runtime(self):
        """
        Get the runtime from the r-* file
        """
        filenames = sorted(["%s/%s"%(self.folder,f) for f in self.run])
        files = [open(filename) for filename in filenames]

        #empty timing
        timing = dict()

        #iterate over all the files
        for f,filename in zip(files,self.run):
            #dictionary for each filename
            this_timing = timing[filename] = dict()
            
            category = "UNKNOWN"
            for line in files[-1]:
                if 'Timing' in line:
                    this_timing[category] = line.split()[-1].split('/')
                if re.search('(\[[0-9.]+\].[A-Z])', line):
                    category = line.strip()

        #close attl the files
        for f in files: f.close()
        self.runtime = timing
        return timing

    def get_data(self,tags):
        """
        Search for a tag in the output files and obtain the data
        """
        data = {}
        for key in self.data.keys():
            if all(tag in key for tag in tags):
                data[key] = dict(list(zip(self.tags[key],np.array(self.data[key]).T)))
        return data

    def print_runtime(self):
        """
        Print the runtime in a string
        """
        timing = self.get_runtime()
        for t in list(timing.items()):
            print(t[0], '\n', t[1], '\n')

    def pack(self,filename=None):
        """
        Pack up all the data in the structure in a json file
        """
        #if no filename is specified we use the same name as the folder
        if not filename: filename = self.folder
        
        #create json dictionary
        jsondata = {"data"     : self.data,
                    "tags"     : self.tags,
                    "runtime"  : self.runtime,
                    "inputfile": self.inputfile,
                    "lattice"  : self.lat,
                    "alat"     : self.alat,
                    "kpts_iku" : self.kpts_iku,
                    "sym_car"  : self.sym_car,
                    "atompos"  : self.atomic_positions,
                    "atomtype" : self.atomic_number}

        #put the data in the file
        filename = '%s.json'%filename
        JsonDumper(jsondata,filename)

    def __str__(self):
        s = ""
        s+= "\nlogs:\n"
        s+= "\n".join(self.logs)+"\n"
        s+= "\nrun:\n"
        s+= "\n".join(self.run)+"\n"
        s+= "\noutput (text):\n"
        s+= "\n".join(self.output)+"\n"
        s+= "\noutput (netcdf):\n"
        s+= "\t\t\n".join(self.netcdf)+"\n"
        s+= "\nlattice:\n"
        s+= "\n".join([("%12.8lf "*3)%tuple(vec) for vec in self.lat])+"\n"
        s+= "\natom positions:\n"
        for an,pos in zip(self.atomic_number,self.atomic_positions):
            s+= "%3d "%an+("%12.8lf "*3)%tuple(pos)+"\n"
        return s

