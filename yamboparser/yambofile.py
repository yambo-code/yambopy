# Copyright (C) 2016 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yamboparser
#
#
import os
import re
import numpy as np

#we try to use netcdf
try:
    from netCDF4 import Dataset
except ImportError:
    _has_netcdf = False
else:
    _has_netcdf = True

class YamboFile():
    """
    This is the Yambo file class.
    It takes as input a filename produced by yambo.
    Can be a netcdf or a text file
    
    List of supported NETCDF files:
        -> ndb.QP

    List of supported text files:
        -> r-*_em?1_*_gw0
        -> o-*.qp
    """
    _output_prefixes = ['o-']
    _report_prefixes = ['r-','r.']
    _log_prefixes    = ['l-','l.']
    _netcdf_prefixes = ['ns','ndb']
    _netcdf_sufixes  = {'QP':'gw'}

    def __init__(self,filename,folder='.'):
        self.filename = filename
        self.folder   = folder
        self.type     = None   
        self.errors   = [] #list of errors
        self.data     = {} #dictionary containing all the important data from the file
 
        if any(prefix in filename for prefix in self._output_prefixes):
            #read lines from file
            f = open("%s/%s"%(folder,filename),'r')
            self.lines = f.readlines()
            f.close()

            #get the line with the title
            try:
                title = self.lines[14]
            except:
                self.errors.append('error reading title')
                return
            if 'GW' in title:
                 self.type = 'output_gw'

        elif any(prefix in filename for prefix in self._report_prefixes):
            self.type = 'report'
        elif any(prefix in filename for prefix in self._log_prefixes):
            self.type = 'log'
        elif any(prefix in filename for prefix in self._netcdf_prefixes):
            for sufix in self._netcdf_sufixes:
                if sufix in filename: 
                    self.type = 'netcdf_%s'%self._netcdf_sufixes[sufix]
                    break

        if self.type is None: self.type = 'unknown'
        
        #parse the file
        self.parse()

    def parse(self):
        """ Parse the file
            Add here things to read log and report files...
        """
        if   self.type == 'netcdf_gw': self.parse_netcdf_gw()
        elif self.type == 'output_gw': self.parse_output()

    def parse_output(self):
        """ Parse an output file from yambo
        """
        #get the tags of the columns
        if self.type == "output_absorption":
            tags = [tag.strip() for tag in re.findall('([ `0-9a-zA-Z\-\/]+)\[[0-9]\]',''.join(self.lines))]
        if self.type == "output_gw":
            tags = [line for line in self.lines if all(tag in line for tag in ['K-point','Band','Eo'])][0]
            tags = tags[2:].strip().split()
        table = np.genfromtxt(self.lines)

        
        self.data = dict(zip(tags,table.T))

    def parse_netcdf_gw(self):
        """ Parse the netcdf gw file
        """
        f = Dataset('%s/%s'%(self.folder,self.filename))
        #quasiparticles table
        qp_table  = f.variables['QP_table'][:].T
        self.data['Kpoint_index'] = qp_table[2]
        self.data['Band'] = qp_table[0]

        #qpoints
        self.data['Kpoint']   = f.variables['QP_kpts'][:].T

        #quasi-particles
        qp = f.variables['QP_E_Eo_Z'][:]
        qp = qp[0]+qp[1]*1j
        self.data['E'], self.data['Eo'], self.data['Z'] = qp.T
        f.close()
        
    def parse_report(self):
        """ Parse the report files.
        Here we check if there were errors running yambo
        """
        pass

    def get_type(self):
        """ Get the type of file        
        """
        return self.type

    def has_errors(self):
        #check if the list is empty
        return not not self.errors

    def get_errors(self):
        """ Check if this is a report file and if it contains errors
        """
        if self.type == 'report':
            return self.errors
        return False

    def get_data(self):
        """ Get the data from this file as a dictionary 
        """
        pass

    def __bool__(self):
        if self.type == 'unknown':
            return False
        else:
            return True
    __nonzero__=__bool__

    def __str__(self):
        return "type: %9s   file: %s/%s"%(self.type, self.folder, self.filename)

