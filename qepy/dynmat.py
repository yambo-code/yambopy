# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
import re
from math import sqrt

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
        s = '&input'
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
