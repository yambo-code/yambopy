from __future__ import print_function
from __future__ import division
# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
from builtins import str
from builtins import zip
from builtins import range
from past.utils import old_div
from builtins import object
import os
import re
from math import sqrt
from numpy import array
from   qepy.auxiliary import *

meVtocm = 8.06573

class DynmatIn(object):
    """
    Generate an manipulate quantum espresso input files for matdyn.x
    """
    def __init__(self):
        self.variable = dict()
        self.qpoints     = [] 

    def write(self,filename):
        f = open(filename,'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        s = '&input\n'
        s += self.stringify_group('',self.variable) #print variable
        if len(self.qpoints) > 0:
          s+=("%d\n"%len(self.qpoints))
          for q in self.qpoints:
            s+=("%12.8lf %12.8lf %12.8lf")%tuple(q[:3])+"\n"
        return s

    def __setitem__(self,key,value):
        self.variable[key] = value

    def __getitem__(self,key):
        return self.variable[key]

    def stringify_group(self, keyword, group):
        if group != {}:
            string='\n'
            for keyword in group:
                string += "%20s = %s\n" % (keyword, group[keyword])
            string += "/\n"
            return string
        else:
            return ''


class Matdyn(object):
    """ Class to read and plot the data from matdyn.modes files 
    """
    _datafile = 'matdyn.modes'
    def __init__(self,natoms,path,folder='.'):
        self.folder   = folder
        self.path     = path
        data_phon     = open("%s/%s"%(folder, self._datafile),'r').readlines()
        self.nmodes   = natoms*3
        self.nqpoints = len(path.get_klist()) 
        self.eigen, self.modes = [], []
        self.qpoints  = []
        for j in range(self.nqpoints):
          frec, v_frec = [], []
          k=2 + j*(self.nmodes*(natoms+1)+5)
          self.qpoints.append(float_from_string(data_phon[k]))
          for i in range(self.nmodes):
            k=4 + j*(self.nmodes*(natoms+1)+5) + i*(natoms+1)
            y = float_from_string(data_phon[k])
            v_mode = []
            for ii in range(1,natoms+1):
              z      = float_from_string(data_phon[k+ii])
              v_atom = array([complex(z[0],z[1]),complex(z[2],z[3]),complex(z[4],z[5])])
              v_mode.append(v_atom)
            v_frec.append(array(v_mode))
            frec.append(y[1])
          self.eigen.append(frec)
          self.modes.append(array(v_frec))

    def plot_eigen(self,path=[]):
        """ plot the phonon frequencies using matplotlib
        """
        import matplotlib.pyplot as plt

        if self.eigen is None:
            print('Error')
            #self.get_eigen()
      
        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            plt.xticks( *list(zip(*path)) )
        plt.ylabel('\\omega (cm$^{-1}$)')

        #plot vertical line
        for point in path:
            x, label = point
            plt.axvline(x)

        #plot bands
        eigen = array(self.eigen)
        for ib in range(self.nmodes):
           plt.plot(range(self.nqpoints),eigen[:,ib], 'r-', lw=2)
        plt.show()
  
    def __str__(self):
        s = ''
        for nq in range(self.nqpoints):
            s+="\n\n q = "+("%12.8lf "*3)%tuple(self.qpoints[nq])+"\n"
            for n in range(self.nmodes):
                s+= "freq (cm-1): %4.3lf\n"%self.eigen[nq][n]
                for na in range(old_div(self.nmodes,3)):
                    xr = self.modes[nq][n][na].real
                    xi = self.modes[nq][n][na].imag
                    s+=("%12.8lf %12.8lfj    "*3)%(xr[0],xi[0],xr[1],xi[1],xr[2],xi[2])+"\n"
        return s
    
    def write_freq_file(self,filename='freq.dat'):
        f = open(filename,'w') 
        for n in range(self.nmodes):
          for nq in range(self.nqpoints):
            f.write("%4.3lf   %4.3lf\n"%(float(nq),self.eigen[nq][n])) 
          f.write("\n") 
        f.close()
