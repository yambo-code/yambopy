# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from subprocess import Popen, PIPE
from yambopy.inputfile import YamboIn
from copy import *
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

class YamboOut():
    """ Class to read yambo output files and pack them in a JSON file

        Input:
        The relative path of the folder where yambo dumped its input files
    """
    _lock = "lock" #name of the lockfile 

    def __init__(self,folder):
        self.folder = folder

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

        self.output = ["%s/%s"%(folder,f) for f in outdir if f[:2] == 'o-']
        self.run    = ["%s/%s"%(folder,f) for f in outdir if f[:2] == 'r-']
        self.logs   = ["%s/LOG/%s"%(folder,f) for f in logdir]
        self.data      = None
        self.runtime   = None
        self.inputfile = None

    def get_data(self):
        """ Search for a tag in the output files and get the data
        """
        files = [open(f) for f in self.output]
        self.data = dict([(filename,np.loadtxt(f)) for filename,f in zip(self.output,files)])
        for f in files: f.close()
        return self.data

    def get_inputfile(self):
        """ Get the input file from the r-* file
        """
        f = open(self.output[-1],'r')
        self.inputfile = []
        for line in f:
            if 'Input file :' in line:
                for line in f:
                    # Note this: to read the input file we just ignore the first 4 characters
                    # of the section after the tag 'Input file:'
                    self.inputfile.append( line[4:] )
        f.close()
        self.inputfile =''.join( self.inputfile )

    def get_runtime(self):
        """ Get the runtime from the r-* file
        """
        files = sorted([open(f) for f in self.run])
        if len(files) > 1: print 'WARNING: more than one file is present, we use the last one (alfabetic order)'
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

    def locked(self):
        """ check if there is a lock present in the folder
        """
        return self._lock in os.listdir(self.folder)

    def put_lock(self):
        """ Put a file to lock the folder.
            This way we don't read anything from this folder again
        """
        f = open(self.folder+'/%s'%self._lock,'w')
        f.close()

    def print_runtime(self):
        """ Print the runtime in a string
        """
        timing = self.get_runtime()
        for t in timing.items():
            print t[0], '\n', t[1], '\n'

    def pack(self,filename=None):
        """ Pack up all the data in the structure in a json file
        """
        if not self.runtime  : self.get_runtime()
        if not self.data     : self.get_data()
        if not self.inputfile: self.get_inputfile()
        if not filename: filename = self.folder

        f = open(filename+'.json','w')
        y = YamboIn()
        json.dump({"data"     : dict(zip(self.data.keys(),[d.tolist() for d in self.data.values()])),
                   "runtime"  : self.runtime,
                   "inputfile": y.read_string(self.inputfile)},
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
        return s


