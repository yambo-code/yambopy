# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from yambopy import *
from yambopy.plot import *
import os
import re
from copy import *
from netCDF4 import Dataset
from yamboparser import *

class YamboOut():
    """ 
    Class to read yambo output files and pack them in a .json file

    **Arguments:**

    ``folder``:      The relative path of the folder where yambo dumped its input files

    ``save_folder``: The path were the SAVE folder is localized 

    """
    _lock = "lock" #name of the lockfile

    def __init__(self,folder,save_folder='.'):

        self.folder = folder

        #check if the save folder is in save_folder if not try folder
        if not os.path.isdir(save_folder+'/SAVE'):
            self.save_folder = folder
        else:
            self.save_folder = save_folder

        #get the output dir
        if os.path.isdir(folder):
            outdir = os.listdir(folder)
        else:
            raise ValueError( "Invalid folder: %s"%folder )
        if os.path.isdir(folder+"/LOG"):
            logdir = os.listdir(folder+"/LOG")
        else:
            logdir = outdir

        tags = ['refl','eel','eps','qp','sf']
        self.output = ["%s"%f for f in outdir if f[:2] == 'o-' and any([tag in f for tag in tags]) and 'xsf' not in f]
        self.run    = ["%s"%f for f in outdir if f[:2] == 'r-']
        self.logs   = ["/LOG/%s"%f for f in logdir]

        # get data from netcdf file ndb.QP (not useful so far)
        netcdf_tags = ['QP']

        #get data from output file
        self.get_runtime()
        self.get_outputfile()
        self.get_inputfile()
        self.get_cell()

        self.netdata = {}  # read netcdf file from YamboFile
        self.nettags = {}
        self.netval  = {}

        # get output name, open netcdf file (only if QP file exists)
        if os.path.exists('%s/ndb.QP' % folder):
          nameout = self.output[0][2:-3]
          self.netdata[nameout] = YamboFile('ndb.QP',folder=folder) 
          self.nettags[nameout] = self.netdata[nameout].data.keys()
          self.netval[nameout]  = self.netdata[nameout].data.values()
          self.set_data_netcdf(nameout)

        # Search of the ndb.QP files. I give the directory of calculations, not the jobname
#        for f in outdir:
#            if os.path.isdir('%s/%s'%(folder,f)) and not 'SAVE' in f:
#                self.netdata[f] = YamboFile('ndb.QP',folder='%s/%s'%(folder,f)) 
#                self.nettags[f] = self.netdata[f].data.keys()
#                self.netval[f]  = self.netdata[f].data.values()
        #fix data from netcdf in suitable format (remove complex type, etc.) 

        # Read data from netcdf file

    def get_cell(self):
        """ 
        Get information about the unit cell (lattice vectors, atom types, positions,
        kpoints and symmetry operations) from the SAVE folder.
        """
        path = self.save_folder+'/SAVE/ns.db1'
        if os.path.isfile(path):
            #read database
            self.nc_db         = Dataset(path)
            self.lat           = self.nc_db.variables['LATTICE_VECTORS'][:].T
            self.alat          = self.nc_db.variables['LATTICE_PARAMETER'][:].T
            self.sym_car       = self.nc_db.variables['SYMMETRY'][:]
            self.kpts_iku      = self.nc_db.variables['K-POINTS'][:].T
            self.apos          = self.nc_db.variables['ATOM_POS'][:,0,:]
            self.atomic_number = self.nc_db.variables['atomic_numbers'][:].T

        else:
            raise ValueError('Could not find ns.db1 in %s from YamboOut'%path)

    def get_outputfile(self):
        """ 
        Get the data from the o-* files
        """
        #open all the o-* files
        files = [open("%s/%s"%(self.folder,f)) for f in self.output]

        self.data = {}
        self.tags = {}
        for filename,f in zip(self.output,files):
            #get the string with the file data
            try:
                string = f.read()
                tags = [tag.strip() for tag in re.findall('#\n#\s+((?:(?:[`0-9a-zA-Z\-\/\|\\(\)\_[\]]+)\s+)+)#\n\s',string)[0].split()]
                f.seek(0)
                self.data[filename] = np.loadtxt(f)
                self.tags[filename] = tags
            except:
                raise ValueError('Error reading file %s'%"%s/%s"%(self.folder,filename))
        #close all the files
        for f in files: f.close()
        return self.data

    def set_data_netcdf(self,nameout):

        self.dictag = {}
        self.dicnet = {}
        self.newtag = []

        val_aux = []

        for word in self.netdata[nameout].data.keys():
            if word == 'Eo':
              self.newtag.append(word)
              val_aux.append( self.netdata[nameout].data[word].real.tolist() ) 
            elif word == 'E':
              self.newtag.append(word)
              val_aux.append( self.netdata[nameout].data[word].real.tolist() ) 
              self.newtag.append('Width[eV]')
              val_aux.append( self.netdata[nameout].data[word].imag.tolist() ) 
            elif word == 'qp_table':
              pass
            elif word == 'E-Eo':
              self.newtag.append(word)
              val_aux.append( self.netdata[nameout].data[word].real.tolist() )
            elif word == 'Band':
              self.newtag.append(word)
              val_aux.append( self.netdata[nameout].data[word].tolist())
            elif word == 'Kpoint_index':
              self.newtag.append(word)
              val_aux.append( self.netdata[nameout].data[word].tolist())
            elif word == 'Kpoint':
              pass
            elif word == 'Z':
              self.newtag.append('Z(Re)')
              val_aux.append( self.netdata[nameout].data[word].real.tolist() )
              self.newtag.append('Z(Im)')
              val_aux.append( self.netdata[nameout].data[word].imag.tolist() )

        # Order elements in a list of variables for each QP states (probably there is a better way)

        n_var = len(val_aux)
        n_qp  = len(val_aux[0])
        aux2 = []
        for i in range(n_qp):
          aux = []
          for j in range(n_var):
            aux.append( val_aux[j][i]  )
          aux2.append(aux)

        # Create a dictionary for tags and another for data

        self.dictag[nameout] = self.newtag
        self.dicnet[nameout] = aux2 

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
        files = sorted([open("%s/%s"%(self.folder,f)) for f in self.run])

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
                data[key] = dict(zip(self.tags[key],np.array(self.data[key]).T))
        return data

    def plot(self,tag,cols=(2,),xlabel=None):
        """
        Search in the output files a certain tag and plot it
        """
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
        """
        Print the runtime in a string
        """
        timing = self.get_runtime()
        for t in timing.items():
            print t[0], '\n', t[1], '\n'

    def pack(self,filename=None):
        """
        Pack up all the data in the structure in a json file
        """
        #if no filename is specified we use the same name as the folder
        if not filename: filename = self.folder
        jsondata = {"data"     : dict(zip(self.data.keys(),[d.tolist() for d in self.data.values()])),
                    "tags"     : self.tags,
                    "runtime"  : self.runtime,
                    "inputfile": self.inputfile,
                    "lattice"  : self.lat,
                    "alat"     : self.alat,
                    "kpts_iku" : self.kpts_iku,
                    "sym_car"  : self.sym_car,
                    "atompos"  : self.apos,
                    "atomtype" : self.atomic_number}
        print (filename)
        filename = '%s.json'%filename
        JsonDumper(jsondata,filename)

    def pack_from_netcdf(self,filename=None):
        """
        Pack up all the data in the structure in a json file
        """
        #if no filename is specified we use the same name as the folder
        if not filename: filename = self.folder
        jsondata = {"data"     : dict(zip(self.dicnet.keys(),self.dicnet.values())),
                    "tags"     : self.dictag,
                    "runtime"  : self.runtime,
                    "inputfile": self.inputfile,
                    "lattice"  : self.lat,
                    "alat"     : self.alat,
                    "kpts_iku" : self.kpts_iku,
                    "sym_car"  : self.sym_car,
                    "atompos"  : self.apos,
                    "atomtype" : self.atomic_number}
        filename = '%s.json'%filename
        JsonDumper(jsondata,filename)

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

