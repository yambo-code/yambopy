# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from subprocess import Popen, PIPE
from yambopy.inputfile import YamboIn
from copy import *
try:
    from netCDF4 import Dataset
    _has_netcdf = True
except ImportError:
    _has_netcdf = False
import os
import json
import numpy as np
import re

#we try to use matplotlib, if not present we won't use it
try:
    from matplotlib import pyplot as plt
except ImportError:
    _has_matplotlib = False
else:
    _has_matplotlib = True

def pack_files_in_folder(folder,save_folder=None):
    """ Helper funciton to look for output files in a
    """
    if not save_folder: save_folder = folder
    #pack the files in .json files
    for dirpath,dirnames,filenames in os.walk(folder):
        #check if there are some output files in the folder
        if ([ f for f in filenames if 'o-' in f ]):
            y = YamboOut(dirpath,save_folder=save_folder)
            y.pack()

class YamboOut():
    """ Class to read yambo output files and pack them in a JSON file

        Input:
        The relative path of the folder where yambo dumped its input files
    """
    _lock = "lock" #name of the lockfile

    def __init__(self,folder,save_folder='./'):
        self.folder = folder
        self.save_folder = save_folder

        #get the output dir
        if os.path.isdir(folder):
            outdir = os.listdir(folder)
        else:
            print "Invalid folder: %s"%folder
            exit(1)
        if os.path.isdir(folder+"/LOG"):
            logdir = os.listdir(folder+"/LOG")
        else:
            logdir = outdir

        self.output = ["%s"%f for f in outdir if f[:2] == 'o-' and ('eel' in f or 'eps' in f or 'qp' in f)]
        self.run    = ["%s"%f for f in outdir if f[:2] == 'r-']
        self.logs   = ["/LOG/%s"%f for f in logdir]
        self.get_runtime()
        self.get_outputfile()
        self.get_inputfile()
        self.get_cell()

    def get_cell(self):
        """ Get information about the unit cell (lattice vectors, atom types, positions,
            kpoints and symmetry operations) from the SAVE folder.
        """
        path = self.save_folder+'/SAVE/ns.db1'
        if os.path.isfile(path) and _has_netcdf:
            #read database
            self.nc_db         = Dataset(path)
            self.lat           = self.nc_db.variables['LATTICE_VECTORS'][:].T
            self.alat          = self.nc_db.variables['LATTICE_PARAMETER'][:].T
            self.sym_car       = self.nc_db.variables['SYMMETRY'][:]
            self.kpts_iku      = self.nc_db.variables['K-POINTS'][:].T
            self.apos          = self.nc_db.variables['ATOM_POS'][:,0,:]
            self.atomic_number = self.nc_db.variables['atomic_numbers'][:].T

        else:
            if not _has_netcdf: print('YamboOut withouth netCDF4 support won\'t retrieve information about the structure')
            print('Could not find ns.db1 in %s'%self.save_folder+'/SAVE')
            self.lat = np.array([])
            self.alat = np.array([])
            self.sym_car = np.array([])
            self.kpts_iku = np.array([])
            self.apos = np.array([])
            self.atomic_number = np.array([])

    def get_outputfile(self):
        """ Get the data from the o-* file
        """
        files = [open("%s/%s"%(self.folder,f)) for f in self.output]
        self.data = dict([(filename,np.loadtxt(f)) for filename,f in zip(self.output,files)])
        for f in files: f.close()
        return self.data

    def get_inputfile(self):
        """ Get the input file from the o-* file
        """
        f = open("%s/%s"%(self.folder,self.output[-1]),'r')
        inputfile = []
        for line in f:
            if 'Input file :' in line:
                for line in f:
                    # Note this: to read the input file we just ignore the first 4 characters
                    # of the section after the tag 'Input file:'
                    inputfile.append( line[4:] )
        f.close()
       
        #use YamboIn to read the input file to a list 
        yi = YamboIn()
        self.inputfile = yi.read_string( ''.join(inputfile) )

    def get_runtime(self):
        """ Get the runtime from the r-* file
        """
        files = sorted([open("%s/%s"%(self.folder,f)) for f in self.run])
        if len(files) > 1:
            print( 'WARNING: more than one r-* file is present in %s'%self.folder )
            print( 'We use the last one (alfabetic order) to get the runtime' )
        timing = dict()
        category = "UNKNOWN"
        for line in files[-1]:
            if 'Timing' in line:
                timing[category] = line.split()[-1].split('/')
            if re.search('(\[[0-9.]+\].[A-Z])', line):
                category = line.strip()
        for f in files: f.close()
        self.runtime = timing
        return timing

    def plot(self,tag,cols=(2,),xlabel=None):
        """ Search in the output files a certain tag and plot it
        """
        if 'self.data' not in locals(): self.get_data()
        for key in self.data.keys():
            if tag in key:
                data = self.data[key]
        for col in cols:
            plt.plot(data[:,0],data[:,col-1],label='col %d'%col)
        plt.title(tag)
        if xlabel: plt.xlabel(xlabel)
        plt.legend()
        plt.show()

    def print_runtime(self):
        """ Print the runtime in a string
        """
        timing = self.get_runtime()
        for t in timing.items():
            print t[0], '\n', t[1], '\n'

    def pack(self,filename=None):
        """ Pack up all the data in the structure in a json file
        """
        #if no filename is specified we use the same name as the folder
        if not filename: filename = self.folder

        f = open('%s.json'%filename,'w')
        json.dump({"data"     : dict(zip(self.data.keys(),[d.tolist() for d in self.data.values()])),
                   "runtime"  : self.runtime,
                   "inputfile": self.inputfile,
                   "lattice"  : self.lat.tolist(),
                   "alat"     : self.alat.tolist(), 
                   "kpts_iku" : self.kpts_iku.tolist(),
                   "sym_car"  : self.sym_car.tolist(),
                   "atompos"  : self.apos.tolist(),
                   "atomtype" : self.atomic_number.tolist()},
                   f,indent=5)
        f.close()

    def __str__(self):
        s = ""
        s+= "\nlogs:\n"
        s+= ("%s\n"*len(self.logs))%tuple(self.logs)
        s+= "\nrun:\n"
        s+= ("%s\n"*len(self.run)%tuple(self.run))
        s+= "\noutput:\n"
        s+= ("%s\n"*len(self.output)%tuple(self.output))
        s+= "\nlattice:\n"
        s+= "\n".join([("%12.8lf "*3)%tuple(vec) for vec in self.lat])+"\n"
        s+= "\natom positions:\n"
        s+= "\n".join(["%3d "%self.atomic_number[n]+("%12.8lf "*3)%tuple(vec) for n,vec in enumerate(self.apos)])+"\n"
        return s
