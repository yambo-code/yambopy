# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
import os
import re
from math import sqrt

__all__ = ['PhIn']

class PhIn(object):
    """
    A class to generate an manipulate quantum espresso input files for ph.x
    """
    def __init__(self):
        self.variable = dict()
        self['tr2_ph']   = 1.0e-12
        self['prefix']   = '\"pw\"'
        self['fildyn']   = '\"dyn\"'
        self['qplot']    = '.false.'
        self['trans']    = '.true.'
        self.qpoints     = [] 

    @classmethod
    def from_qpoints(cls,qpoints):
        """
        Create a ph.x input from qpoints
        If qpoints is a list of lists write those qpoints explicit
        if qpoints is a list of indexes set the nq1,nq2 and nq3 indexes
        """
        instance = cls()
        if isinstance(qpoints[0],list): instance.qpoints = qpoints
        else: instance.set_nq(*qpoints)
        return instance

    @property
    def prefix(self):
        return self['prefix'].replace("'",'')

    @prefix.setter
    def prefix(self,value):
        self['prefix'] = "'%s'"%value.replace("'",'')

    @property
    def fildyn(self):
        return self['fildyn'].replace("'",'')

    @fildyn.setter
    def fildyn(self,value):
        self['fildyn'] = "'%s'"%value.replace("'",'')

    def write(self,filename):
        with open(filename,'w') as f:
            f.write(str(self))

    def set_vars(**kwargs):
        for var,value in kwargs.items():
            self[var] = value

    def set_nq(self,nq1,nq2,nq3):
        self['nq1'] = nq1
        self['nq2'] = nq2
        self['nq3'] = nq3
        self['ldisp'] = ".true."

    def __str__(self):
        lines = []; app = lines.append
        app("\n phonons \n&inputph")
        app(self.stringify_group('',self.variable)) #print variable
        if 'true' in self['qplot'].lower(): 
            app("%d"%len(self.qpoints))
            for q in self.qpoints:
                app("%12.8lf %12.8lf %12.8lf %d"%tuple(q))
        return "\n".join(lines)+"\n"

    def __setitem__(self,key,value):
        self.variable[key] = value

    def __getitem__(self,key):
        return self.variable[key]

    def stringify_group(self, keyword, group):
        if group != {}:
            lines = []; app = lines.append
            for keyword in group:
                app("%10s = %s" % (keyword, group[keyword]))
            app("/")
            return "\n".join(lines)
        else:
            return ''
