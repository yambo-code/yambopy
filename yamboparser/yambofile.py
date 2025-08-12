# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yamboparser
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

def if_has_netcdf(f):
    if _has_netcdf:
        return f
    
class YamboFile(object):
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
    _log_prefixes    = ['l-','l.','l_']
    _netcdf_prefixes = ['ns','ndb']
    _netcdf_sufixes  = {'QP':'gw','HF_and_locXC':'hf'}
    _outputs_type = {'output_abs':'eps','output_loss':'eel','output_alpha':'alpha', 'output_jdos':'jdos'}

    def __init__(self,filename,folder='.',**parse_kwargs):
        self.filename = filename
        self.folder = folder
        self.errors = [] #list of errors
        self.warnings = [] #list of warnings
        self.memstats = [] #list of memory allocation statistics
        self.data = {} #dictionary containing all the important data from the file
        self.kpoints = {}
        self.timing = []

        #get the type of file
        self.type = YamboFile.get_filetype(filename,folder)

        #if needed read the lines
        if self.type in ['output_gw', 'log', 'report', 'output_abs', 'output_loss', 'output_alpha', 'output_jdos',]:
            #read lines from file
            with open(os.path.join(folder,filename),'r') as f:
                self.lines = f.readlines()
            
        #parse the file
        self.parse(**parse_kwargs)
    
    @staticmethod
    def get_filetype(filename,folder):
        """
        Get the type of file
        """
        type = 'unknown'
        basename = os.path.basename(filename)
        if any(basename.startswith(prefix) for prefix in YamboFile._output_prefixes):
            #read lines from file
            with open(os.path.join(folder,basename),'r') as f:
                lines = f.readlines()

            #get the line with the title
            title = lines[14]

            if 'GW' in title or '.qp' in filename:
                 type = 'output_gw'
            
            else:
                for key, val in list(zip(YamboFile._outputs_type.keys(),YamboFile._outputs_type.values())):
                    if val in filename: type = key

        elif any(basename.startswith(prefix) for prefix in YamboFile._report_prefixes):
            type = 'report'
        elif any(basename.startswith(prefix) for prefix in YamboFile._log_prefixes):
            type = 'log'
        elif any(basename.startswith(prefix) for prefix in YamboFile._netcdf_prefixes):
            for sufix in YamboFile._netcdf_sufixes:
                if basename.endswith(sufix):
                    type = 'netcdf_%s'%YamboFile._netcdf_sufixes[sufix]
                    break

        return type

    def parse(self,**parse_kwargs):
        """ Parse the file
            Add here things to read log and report files...
        """
        if   self.type == 'netcdf_gw': self.parse_netcdf_gw(**parse_kwargs)
        elif self.type == 'netcdf_hf': self.parse_netcdf_hf(**parse_kwargs)
        elif self.type in ['output_gw', 'output_abs', 'output_loss', 'output_alpha', 'output_jdos']: self.parse_output(**parse_kwargs)
        elif self.type == 'log': self.parse_log(**parse_kwargs)
        elif self.type == 'report'  : self.parse_report(**parse_kwargs)

    def parse_output(self,**parse_kwargs):
        """ Parse an output file from yambo,
        """
        zip_tags = parse_kwargs.get('zip_tags',False) #flag--default behavior is to do nothing
        #get the tags of the columns
        if self.type in YamboFile._outputs_type.keys():  #== "output_absorption":
            pattern = '([ `0-9a-zA-Z\-\/]+)\[[0-9]\]'
            tags = [tag.strip() for tag in re.findall(pattern,''.join(self.lines))]
            if len(tags) < 2 and self.type in ["output_loss", "output_abs", "output_alpha"]: # temporary fix for IP case: E[1] [eV]          Im(eps)            Re(eps)
                lines_with_matches = [line for line in self.lines if re.search(pattern, line)] 
                tags = [tag.strip() for tag in lines_with_matches[0].replace("#","").replace("\n","").replace("[eV]","").split()]
        if self.type == "output_gw":
            tags = [line.replace('(meV)','').replace('Sc(Eo)','Sc|Eo') for line in self.lines if all(tag in line for tag in ['K-point','Band','Eo'])][0]
            tags = tags[2:].strip().split()
        table = np.loadtxt(self.lines)
        _kdata ={}
        k_index =[ str(int(i)) for i in table[:,0]] # first column  has kpoints
        for ind in range(len(k_index)):
            for itag in range(len(tags)):
                 if k_index[ind] not in list(_kdata.keys()):
                     _kdata[k_index[ind]] = {}
                 try:
                     _kdata[k_index[ind]][tags[itag]].append(table[ind,itag]) #errors when you have multiple qp? IndexError: index 5 is out of bounds for axis 1 with size 5
                 except KeyError:
                     _kdata[k_index[ind]][tags[itag]]  = [ table[ind,itag] ]

        self.data = _kdata
        if (zip_tags): #combines tags such that keys refer to the columns in data file
            self.data = dict(zip(tags,table.T))

    @if_has_netcdf
    def parse_netcdf_gw(self,**parse_kwargs):
        """ Parse the netcdf gw file
        """
        data = {}

        filename = '%s/%s'%(self.folder,self.filename)
        with Dataset(filename) as f:

            #quasiparticles table
            qp_table  = f.variables['QP_table'][:]
            data['Kpoint_index'] = qp_table[2]
            data['Band'] = qp_table[0]
            #print(qp_table.shape)
            if qp_table.shape[0] == 4: # spin polarized
                data['Spin_pol'] = qp_table[3]
            data['qp_table'] = qp_table[:]  # ib, ik, ,(isp if spin polarized)
            #qpoints
            data['Kpoint']   = f.variables['QP_kpts'][:].T

            #quasi-particles
            if 'QP_E_Eo_Z' in f.variables:
                #old format
                qp = f.variables['QP_E_Eo_Z'][:]
                qp = qp[0]+qp[1]*1j
                data['E'],  data['Eo'], data['Z'] = qp.T
                data['E-Eo'] = data['E']  -  data['Eo']
                self.data=data
            else:
                #new format
                E  = f.variables['QP_E'][:]
                data['E'] = E[:,0] + E[:,1]*1j
                Eo = f.variables['QP_Eo'][:]
                data['Eo']= np.array(Eo,dtype=data['E'].dtype)
                Z  = f.variables['QP_Z'][:]
                data['Z'] = Z[:,0] + Z[:,1]*1j
                data['E-Eo'] = data['E']  -  data['Eo']
                self.data=data

    @if_has_netcdf
    def parse_netcdf_hf(self,**parse_kwargs):
        """ Parse the netcdf hf file (ndb.HF_and_locXC)
        """
        data = {}

        f = Dataset(os.path.join(self.folder,self.filename))

        pars = f.variables['PARS'][:]

        data['nkpoints'] = int(pars[0])
        data['nbands'] = int(pars[1])
        #data['QP_table'] = f.variables['QP_table'][:]

        #old format
        if 'Sx_Vxc' in f.variables:
            hf = f.variables['Sx_Vxc'][:]
            if hf.shape[0]%8 ==0 :
                qp =  hf.reshape(-1,8)
                ib, ibp, ik, isp, rsx, isx, revx, imvx = qp.T
            else:
                qp =  hf.reshape(-1,7)
                ib, ibp, ik, rsx, isx, revx, imvx = qp.T
            data['Sx'] = rsx + isx*1j
            data['Vxc'] = revx + imvx*1j
        #new format
        else:
            Sx  = f.variables['Sx'][:]
            data['Sx'] = Sx[:,0] + Sx[:,1]*1j
            Vxc = f.variables['Vxc'][:]
            data['Vxc'] = Vxc[:,0] + Vxc[:,1]*1j
            
        self.data=data
        f.close()

    def parse_report(self,**parse_kwargs):
        """ Parse the report files.
            produces output of this nature:
            { k-index1  : { 'dft_enrgy':[...], 'qp_energy':[...] },
              k-index2  :{...}
            }
            k-index is the kpoint at which the yambo calculation was
            done.
        """
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

    def parse_log(self,**parse_kwargs):
        """ Get ERRORS and WARNINGS from  l-*  file, useful for debugging
        """
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

