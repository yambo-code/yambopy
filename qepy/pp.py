# Author = Riccardo Reho r.reho@uu.nl
#
# This file is part of yambopy
#
# Interface with pp.x, post-processing tool of QE

import os
import re
import shutil
import numpy as np
from qepy import qepyenv

class PPIn():
    """
    Class to generate an manipulate Quantum Espresso post-processing pp.x input files
    This class is not meant to be comprehensive but a lightweight version 
    capable of basic input/ouput of PW files.

    Examples of use:

    To write a local file from scratch with name "pp.in"

        .. code-block :: python

            pp = PwIn.from_file('prefix.pp')
            print pp

    For now, this has been tested only for the plotting of wavefunctions |psi**2| from QE, as
    in the followin example.

        .. code-block :: python

            from yambopy import *
            prefix = 'prefix'
            kpoints =np.arange(1,2,1)
            kbands = np.arange(456,457,1)
            nf = (kbands[-1]-kbands[0]+1)*len(kpoints) #+1 due to the usual python counting
            print ("Post-processing pp.x\n nfile =", nf,"\n")

            files = []
            weights=np.ones(nf)
            for i_k, k in enumerate(kpoints):
                for i_b,b in enumerate(kbands):
                    files.append(f"'./wfc-datfiles/{prefix}_K{str(k).zfill(3)}_B{str(b).zfill(3)}'")
            # Initialize PPIn class
            pp = PPIn()
            # In this specific case, the .dat files are already available from previous run, therefore I only need to plot
            pp.set_nfile(nf)
            pp.set_files(files)
            pp.set_weigths(weights)
            pp.set_fileout(f"'plot-B{kbands[-1]}-B{kbands[0]}-K{kpoints[-1]}-K{kpoints[0]}.xsf'")
            pp.write(f'plot-B{kbands[0]}-B{kbands[-1]}-K{kpoints[0]}-K{kpoints[-1]}.in',plotonly=True)

    """

    _pp = "pp.x"

    def __init__(self):
        """ TODO: specify the required parameters """
        #dictionaries
        self.inputpp = dict()
        self.plot = dict(nfile = 1, output_format=5, fileout = "'plot.xsf'", iflag=3)

    def set_prefix(self,prefix):
        self.inputpp['prefix'] = prefix
    
    def set_outdir(self,outdir):
        self.inputpp['outdir'] = outdir

    def set_filplot(self,filplot):
        self.inputpp['filplot'] = filplot

    def set_plotnum(self,plotnum):
        self.inputpp['plotnum'] = plotnum

    def set_kpoints(self,kpoints):
        "k-point should be an array with a range of values you want to compute the quantities or a single k-point"
        if (len(kpoints)==1):
            if (isinstance(kpoints,int)):
                self.inputpp['kpoint(1)'] = kpoints
            if (isinstance(kpoints,list)):
                self.inputpp['kpoints(1)'] = kpoints[0]
        else: 
            if(len(kpoints) >1):
                self.inputpp['kpoints(1)'] = kpoints[0]
                self.inputpp['kpoints(2)'] = kpoints[1]

    def set_kbands(self,kbands):
        "k-point should be an array with a range of values you want to compute the quantities or a single k-point"
        if (len(kbands)==1):
            if (isinstance(kbands,int)):
                self.inputpp['kbands(1)'] = kbands
            if (isinstance(kbands,list)):
                self.inputpp['kbands(1)'] = kbands[0]
        else:
            if(len(kbands) >1):
                self.inputpp['kbands(1)'] = kbands[0]
                self.inputpp['kbands(2)'] = kbands[1]

    def set_nfile(self,nfile):
        self.plot['nfile'] = nfile

    def set_files(self,files):
        for i_f,f in enumerate(files):
            self.plot[f'filepp({i_f+1})'] = f #+1 for python-Fortran counting
    
    def set_weigths(self,weights):
        for i_w,w in enumerate(weights):
            self.plot[f'weight({i_w+1})'] = w #+1 for python-Fortran counting

    def set_outputformat(self,outputformat):
        self.plot['output_format'] = outputformat
    
    def set_fileout(self,fileout):
        self.plot['fileout'] = fileout

    def set_iflag(self,iflag):
       self.plot['iflag'] = iflag
       
    def stringify_group(self, keyword, group):
        if group != {}:
            string='&%s\n' % keyword
            for keyword in sorted(group): # Py2/3 discrepancy in keyword order
                string += "    %s = %s\n" % (keyword, group[keyword]) #4 space are mandatory at the beg of the line for pp.x syntax
            string += "/"
            return string
        else:
            return ''
    
    def write(self,filename,plotonly=False):
        """write the file to disk """
        with open(filename,'w') as f:
           if (plotonly):
              """an empyt block for inputpp is mandatory in QE if you only want to plot results from data"""
              f.write('&inputpp \n/\n')
           f.write(str(self))               
    def get_string(self):
        """
        Output the file in the form of a string
        """
        lines = []; app = lines.append
        app( self.stringify_group("inputpp",self.inputpp) ) #print inputpp
        app( self.stringify_group("plot",self.plot) ) #print plot
    
        return "\n".join(lines)

    def __str__(self):
        return self.get_string()    
