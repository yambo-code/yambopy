# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
import re
from math import sqrt

class PhIn():
    """ A class to generate an manipulate quantum espresso input files for ph.x
    """
    def __init__(self):
        self.variable = dict()
        self['tr2_ph']   = 1.0e-12
        self['prefix']   = '\"pwscf\"'
        self['fildyn']   = '\"matdyn\"'
        self['fildvscf'] = '\"\"'
        self['qplot']    = '.false.'
        self['trans']    = '.true.'
        self.qpoints     = [] 

    def write(self,filename):
        f = open(filename,'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        s = '&inputph'
        s += self.stringify_group('',self.variable) #print variable
        if 'true' in self['qplot'].lower(): 
            s += "%d\n"%len(self.qpoints)
            for q in self.qpoints:
                s+=("%12.8lf %12.8lf %12.8lf %d")%tuple(q)+"\n"
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
