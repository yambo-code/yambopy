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
        self.warnings   = [] #list of warnings
        self.memstats   = [] #list of memory allocation statistics
        self.data     = {} #dictionary containing all the important data from the file
        self.kpoints = {}
        self.timing = []
 
        if any(filename.startswith(prefix) for prefix in self._output_prefixes):
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

        elif any(filename.startswith(prefix) for prefix in self._report_prefixes):
            self.type = 'report'
        elif any(filename.startswith(prefix) for prefix in self._log_prefixes):
            self.type = 'log'
        elif any(filename.startswith(prefix) for prefix in self._netcdf_prefixes) and _has_netcdf:
            for sufix in self._netcdf_sufixes:
                if filename.endswith(sufix): 
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
        elif self.type == 'log': self.parse_log()
        elif self.type == 'report'  : self.parse_report()

    def parse_output(self):
        """ Parse an output file from yambo,
        """
        #get the tags of the columns
        if self.type == "output_absorption":
            tags = [tag.strip() for tag in re.findall('([ `0-9a-zA-Z\-\/]+)\[[0-9]\]',''.join(self.lines))]
        if self.type == "output_gw":
            tags = [line for line in self.lines if all(tag in line for tag in ['K-point','Band','Eo'])][0]
            tags = tags[2:].strip().split()
        table = np.genfromtxt(self.lines)
        _kdata ={}
        k_index =[ str(int(i)) for i in table[:,0]] # first column  has kpoints
        for ind in range(len(k_index)):
            for itag in range(len(tags)):
                 if k_index[ind] not in _kdata.keys():
                     _kdata[k_index[ind]] = {}
                 try:
                     _kdata[k_index[ind]][tags[itag]].append(table[ind,itag])
                 except KeyError:
                     _kdata[k_index[ind]][tags[itag]]  = [ table[ind,itag] ]

        self.data = _kdata 
        #self.data = dict(zip(tags,table.T))

    def parse_netcdf_gw(self):
        """ Parse the netcdf gw file
        """
        if _has_netcdf:

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
            produces output of this nature:
            { k-index1  : { 'dft_enrgy':[...], 'qp_energy':[...] },
              k-index2  :{...}
            }
            k-index is the kpoint at which the yambo calculation was
            done.
        """
        if not hasattr(self, 'lines'):
            with open('%s/%s'%(self.folder,self.filename)) as fl:
                self.lines = fl.readlines()
        # start with check for  failure due to error:
        err = re.compile('^\s+?\[ERROR\]\s+?(.*)$')
        kpoints = re.compile('^  [A-X*]+\sK\s\[([0-9]+)\]\s[:](?:\s+)?([0-9.E-]+\s+[0-9.E-]+\s+[0-9.E-]+)\s[A-Za-z()\s*.]+[0-9]+[A-Za-z()\s*.]+([0-9.]+)')
        memory = re.compile('^\s+?<([0-9a-z-]+)> ([A-Z0-9]+)[:] \[M  ([0-9.]+) Gb\]? ([a-zA-Z0-9\s.()\[\]]+)?')
        timing = re.compile('\s+?[A-Za-z]+iming\s+?[A-Za-z/\[\]]+[:]\s+?([a-z0-9-]+)[/]([a-z0-9-]+)[/]([a-z0-9-]+)')
        self.memstats.extend([ line for line in self.lines if memory.match(line)])
        for line in self.lines:
            if err.match(line):
                if 'STOP' in err.match(line).groups()[0]:
                    # stop parsing, this is a failed calc.
                    self.errors.append(err.match(line).groups()[0])
                    return 
            if timing.match(line):
                self.timing.append(timing.match(line).groups()[0] )
            if kpoints.match(line):
                kindx, kpt, wgt = kpoints.match(line).groups()
                self.kpoints[str(int(kindx))] =  [ float(i.strip()) for i in kpt.split()]
                    
        full_lines = '\n'.join(self.lines)
        qp_regx = re.compile('(^\s+?QP\s\[eV\]\s@\sK\s\[\d+\][a-z0-9E:()\s.-]+)(.*?)(?=^$)',re.M|re.DOTALL)
        kp_regex = re.compile('^\s+?QP\s\[eV\]\s@\sK\s\[(\d+)\][a-z0-9E:()\s.-]+$')
        spliter = re.compile('^(B[=]\d+\sEo[=]\s+?[E0-9.-]+\sE[=]\s+?[E0-9.-]+\sE[-]Eo[=]\s+?[E0-9.-]+\sRe[(]Z[)][=]\s+?[E0-9.-]+\sIm[(]Z[)][=]\s?[E0-9.-]+\snlXC[=]\s+?[E0-9.-]+\slXC[=]\s+?[E0-9.-]+\sSo[=]\s+?[E0-9.-]+)')
        extract = re.compile('B[=](\d+)\sEo[=](?:\s+)?([E0-9.-]+)\sE[=](?:\s+)?([E0-9.-]+)\sE[-]Eo[=](?:\s+)?([E0-9.-]+)\sRe[(]Z[)][=](?:\s+)?([E0-9.-]+)\sIm[(]Z[)][=](?:\s+)?[E0-9.-]+\snlXC[=](?:\s+)?([E0-9.-]+)\slXC[=](?:\s+)?([E0-9.-]+)\sSo[=](?:\s+)?([E0-9.-]+)')
        qp_lines = qp_regx.findall(full_lines)
        qp_results ={}
        for each in qp_lines: # first group of qp data, shares k-point index
            kp_index = None
            kp_results={'bindex':[],'dft_energy':[],'qp_energy':[],'qp_correction':[],
                        'z_factor':[],'non_local_xc':[],'local_xc':[],'selfenergy_c':[]}
            for line in each: # different band indexes =B
                if kp_regex.match(line):
                    kp_index = str(int(kp_regex.match(line).groups()[0]))
                else: #  data line B=x Eo = y ..
                    data_lines = [ i for i in spliter.split(line) if i.strip()]
                    for qp_data in data_lines:
                        bindex, dft_energy, qp_energy, qp_correction, z_factor, \
                        non_local_xc, local_xc, selfenergy_c = [float (i) for i in extract.match(qp_data).groups()]
                        kp_results['bindex'].append(bindex) 
                        kp_results['dft_energy'].append(dft_energy) 
                        kp_results['qp_energy'].append(qp_energy) 
                        kp_results['qp_correction'].append(qp_correction) 
                        kp_results['z_factor'].append(z_factor) 
                        kp_results['non_local_xc'].append(non_local_xc) 
                        kp_results['local_xc'].append(local_xc) 
                        kp_results['selfenergy_c'].append(selfenergy_c) 
            qp_results[kp_index] = kp_results
        self.data = qp_results 

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

    def parse_log(self):
        """ Get ERRORS and WARNINGS from  l-*  file, useful for debugging
        """
        if not hasattr(self, 'lines'):
            with open('%s/%s'%(self.folder,self.filename)) as fl:
                self.lines = fl.readlines()
        warning = re.compile('^\s+?<([0-9a-z-]+)> ([A-Z0-9]+)[:] \[(WARNING)\]? ([a-zA-Z0-9\s.()\[\]]+)?')
        error = re.compile('^\s+?<([0-9a-z-]+)> ([A-Z0-9]+)[:] \[(ERROR)\]? ([a-zA-Z0-9\s.()\[\]]+)?')
        self.warnings.extend([ line for line in self.lines if warning.match(line)])
        self.errors.extend([ line for line in self.lines if error.match(line)])
        

    def __bool__(self):
        if self.type == 'unknown':
            return False
        else:
            return True
    __nonzero__=__bool__

    def __str__(self):
        return "type: %9s   file: %s/%s"%(self.type, self.folder, self.filename)

