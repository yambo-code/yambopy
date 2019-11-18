# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function
import os

__all__ = ['ProjwfcIn']

class ProjwfcIn():
    """
    A class to generate and manipulate projwfc input files.
    """    
    _projwfc = 'projwfc.x'

    def __init__(self,prefix,DeltaE = 0.05,folder='.'):
        self.prefix = prefix
        self.DeltaE = DeltaE
        
    def run(self,filename=None,procs=1,folder='.'):
        """ 
        Function is used to run this job locally
        """
        filename = self.getfilename(filename)
        os.system('mkdir -p %s'%folder)
        self.write("%s/%s"%(folder,filename))
        print("Run projwfc...")
        if procs == 1:
            os.system('cd %s; %s -inp %s > projwfc.log' % (folder,self._projwfc,filename))
        else:
            os.system('cd %s; mpirun -np %d %s -inp %s > projwfc.log' % (folder,procs,self._projwfc,filename))
        print("done!")

    def getfilename(self,filename):
        if filename is None:
            filename = "%s.projwfc"%self.prefix
        return filename

    def write(self,filename=None,folder='.'):
        """
        Write the projwfc file
        """
        filename = self.getfilename(filename)
        os.system('mkdir -p %s'%folder)
        f = open("%s/%s"%(folder,filename),'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        """
        Output the file in the from of a string
        """
        string = '\n'
        string += "&projwfc\n"
        string += "prefix = '%s'\n"%self.prefix
        string += "DeltaE = %lf\n"%self.DeltaE
        string += "/\n"
        return string
