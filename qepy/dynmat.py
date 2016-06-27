# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
import re
from math import sqrt
from numpy import array
from   qepy.auxiliary import *

meVtocm = 8.06573

class DynmatIn():
    """ A class to generate an manipulate quantum espresso input files for matdyn.x
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


class Matdyn():
    """ Class to read data from matdyn.modes files 
    """
    _datafile = 'matdyn.modes'
    def __init__(self,natoms,nqpoints,path='.'):
        self.path     = path
        data_phon     = open("%s/%s"%(path, self._datafile),'r').readlines()
        self.nmodes   = natoms*3
        self.nqpoints = nqpoints
        self.eigen, self.modes = [], []
        for j in xrange(nqpoints):
          frec, v_frec = [], []
          for i in xrange(self.nmodes):
            k=4 + j*(self.nmodes*(natoms+1)+5) + i*(natoms+1)
            y = float_from_string(data_phon[k])
            v_mode = []
            for ii in xrange(1,natoms+1):
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
###################################################################################
# Font selection and borders
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
###################################################################################

        if self.eigen is None:
            print('Error')
            #self.get_eigen()
      
        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            plt.xticks( *zip(*path) )
        plt.ylabel('\omega (cm$^{-1}$)')

        #plot vertical line
        for point in path:
            x, label = point
            plt.axvline(x)

        #plot bands
        eigen = array(self.eigen)
        for ib in range(self.nmodes):
           plt.plot(xrange(self.nqpoints),eigen[:,ib], 'r-', lw=2)
        plt.show()
  

class Matdyn():
    """ Class to read and plot the data from matdyn.modes files 
    """
    _datafile = 'matdyn.modes'
    def __init__(self,natoms,nqpoints,path='.'):
        self.path     = path
        data_phon     = open("%s/%s"%(path, self._datafile),'r').readlines()
        self.nmodes   = natoms*3
        self.nqpoints = nqpoints
        self.eigen, self.modes = [], []
        for j in xrange(nqpoints):
          frec, v_frec = [], []
          for i in xrange(self.nmodes):
            k=4 + j*(self.nmodes*(natoms+1)+5) + i*(natoms+1)
            y = float_from_string(data_phon[k])
            v_mode = []
            for ii in xrange(1,natoms+1):
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
###################################################################################
# Font selection and borders
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
###################################################################################

        if self.eigen is None:
            print('Error')
            #self.get_eigen()
      
        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            plt.xticks( *zip(*path) )
        plt.ylabel('\omega (cm$^{-1}$)')

        #plot vertical line
        for point in path:
            x, label = point
            plt.axvline(x)

        #plot bands
        eigen = array(self.eigen)
        for ib in range(self.nmodes):
           plt.plot(xrange(self.nqpoints),eigen[:,ib], 'r-', lw=2)
        plt.show()
  
