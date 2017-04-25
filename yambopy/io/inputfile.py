from __future__ import print_function
from past.builtins import basestring
# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
from builtins import str
from builtins import map
from builtins import object
from subprocess import Popen, PIPE
import os
import json
from time import sleep
import re

class YamboIn(object):
    """
    Class to read, write, create and manipulate yambo input files with python.

    Examples of use:

    Initialize an empty input file:

        .. code-block:: python

            y = YamboIn(filename='somefile.in')
            print y

    Call yambo to initialize the input file with variables according to the runlevel,
    parse the input file and store the variables:

        .. code-block:: python

            y = YamboIn('yambo -o c',folder='ip')
            print y

    If the argument ``args`` was used then the filename should be left as ``yambo.in`` because that's the default input filename that yambo will create.

    Call ypp to initialize the input file:

        .. code-block:: python

            y = YamboIn('yyp -e w'args=,filename='ypp.in')
            print y

    **Arguments:**

        ``args``:     if specified yambopy will run yambo, read the generated input file and initialize the class with those variables.

        ``folder``:   the folder where the SAVE directory is located

        ``vim``:      if yambo is compiled using vim as an editor this variable should be set as True because then `yambopy` will close vim.
        In newer versions an argument for yambo '-Q' tells it to not call vim

        ``filename``: the name of the input file to be read

    """
    #Regular expressions
    _variaexp   = '([A-Za-z\_0-9]+(?:\_[A-Za-z]+)?)' #variables names
    _numexp     = '([+-]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)' #number
    _spacexp    = '(?:[ \t]+)?' #space
    _stringexp  = '["\']([a-zA-Z0-9_ ]+?)["\']' #string
    _arrayexp   = '%'+_spacexp+_variaexp+'\s+(?:\#.+)?((?:(?:\s|\.|[+-]?\d)+?\|)+)\s+([a-zA-Z]+)?' #arrays
    _complexexp = '\('+_spacexp+_numexp+_spacexp+','+_spacexp+_numexp+_spacexp+'\)' #complex numbers
    _runexp     = '([a-zA-Z0-9_]+)' #runlevels
    # list of available runlevels to be stored in the arguments array
    _runlevels  = ['rim_cut','chi','em1s','bse','optics','bsk','bss',
                   'em1d','gw0','HF_and_locXC','setup','ppa','cohsex','life',
                   'collisions','negf','el_ph_scatt','el_el_scatt','excitons','wavefunction','fixsyms',
                   'QPDBs', 'QPDB_merge','RealTime','RT_X','RToccDos','RToccBnd','RToccEner',
                   'RToccTime','RTlifeBnd','amplitude','bzgrids','Random_Grid','gkkp','el_ph_corr','WRbsWF','Select_energy', 'RTDBs','photolum','kpts_map',
                   'RTtime','RToccupations','RTfitbands']

    def __init__(self,args='',folder='.',filename='yambo.in'):
        """
        Initalize the class
        """
        self.folder = folder
        self.yamboargs = args

        #the type of the variables is determined from the type of variable in this dictionary
        self.variables = {} #here we will store the values of the variables
        self.arguments = [] #here we will store the arguments

        # if we initalize the class with arguments we call yambo to generate the input file
        if args != '':
            workdir = os.getcwd()
            os.chdir(folder)
            os.system('rm -f %s'%filename)
            yambo = Popen(args+' -Q', stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True)
            print('YAMBOPY: yambo command: %s'%(args+' -Q'))
            yambo.wait()
            os.chdir(workdir)
            self.read_file(filename="%s/%s"%(folder,filename))
        else:
            if filename:
                self.read_file(filename="%s/%s"%(folder,filename))

    def __getitem__(self,key):
        """ Get the value of a variable in the input file
        """
        return self.variables[key]

    def __setitem__(self,key,value):
        """ Set the value of a variable in the input file
        """
        #if the units are not specified, add them
        if isinstance(value,list):
            if not any( [isinstance(v,basestring) for v in value] ):
                value = [value,'']
        if isinstance(value,(int,float,complex)):
            value = [value,'']
        self.variables[key] = value

    def __delitem__(self,key):
        """ Remove a keyword from the dicitonary
        """
        del self.variables[key]

    def read_file(self,filename='yambo.in'):
        """ Read the variables from a file
        """
        try:
            yambofile = open(filename,"r")
        except IOError:
            print('Could not read the file %s'%filename)
            print('Something is wrong, yambo did not create the input file. Or the file you are trying to read does not exist')
            print('command: %s'%self.yamboargs)
            print('folder:  %s/'%self.folder)
            exit()
        inputfile = self.read_string(yambofile.read())
        yambofile.close()

    def add_dict(self,variables):
        """
        Add a dictionary containing variables to the current inputfile
        """
        self.variables.update(variables)

    def read_string(self,inputfile):
        """
        Read the input variables from a string
        """
        var_real     = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp +
                                  self._numexp + self._spacexp + '([A-Za-z]+)?',inputfile)
        var_string   = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp + self._stringexp, inputfile)
        var_array    = re.findall(self._arrayexp,inputfile)
        var_complex  = re.findall(self._variaexp + self._spacexp + '='+ self._spacexp + self._complexexp + self._spacexp + '([A-Za-z]+)?', inputfile)
        var_runlevel = re.findall(self._runexp + self._spacexp, inputfile)

        def clean(a):
            """
            clean the variables according to the type of data
            """
            a = a.strip()
            if a.replace('.','',1).isdigit():
                if "." in a: return float(a)
                else:        return int(a)
            return a

        # Determination of the arguments
        for key in self._runlevels:
            if key in var_runlevel:
                self.arguments.append(key)

        #float variables
        for var in var_real:
            name, value, unit = var
            self[name] = [float(value),unit]

        #string variables
        for var in var_string:
            name, string = var
            self[name] = string

        #complex variables
        for var in var_complex:
            name, real, imag, unit = var
            self[name] = [complex(float(real),float(imag)),unit]

        #array variables
        for var in var_array:
            name, array, unit = var
            array = [clean(val) for val in array.split('|')[:-1]]
            self[name] = [array,unit]

        return {"arguments": self.arguments, "variables": self.variables}

    def optimize(self,conv,variables=('all',),run=lambda x: None,ref_run=True):
        """
        Function to to make multiple runs of yambo to converge calculation parameters

        Arguments:

            A dictionary conv that has all the variables to be optimized
            A list fo the name of the variables in the dicitonary that are to be optimized
            A function run that takes as input the name of the inputfile (used to run yambo)
            A boolean ref_run that can disable the submitting of the reference run (see scripts/analyse_gw.py)

        """
        name_files = []

        #check which variables to optimize
        if 'all' in variables:
            variables = list(conv.keys())

        #save all the variables
        backup = {}
        for var in variables:
            backup[var] = self[var]

        #add units to all the variables (to be improved!)
        for key,value in list(conv.items()):
            if not isinstance(value[-1],basestring) and type(value[0]) == list:
                conv[key] = [value,'']

        #make a first run with all the first elements
        reference = {}
        for key,value in list(conv.items()):
            values, unit = value
            reference[key] = [values[0],unit]
            self[key] = [values[0],unit]
        #write the file and run
        if ref_run==True:
            self.write( "%s/reference.in"%(self.folder) )
            run('reference.in')
        else:
            print('Reference run disabled.')

        #converge one by one
        for key in [var for var in list(conv.keys()) if var in variables]:
            values, unit = conv[key]
            #put back the original values of the variables
            for var in variables:
                self[var] = reference[var]
            #change the other variables
            if isinstance(values[0],basestring):
                for string in values[1:]:
                    filename = "%s_%s"%(key,string)
                    self[key] = string
                    self.write( self.folder + filename )
                    run(filename+".in")
                continue
            if type(values[0])==float:
                for val in values[1:]:
                    filename = "%s_%12.8lf"%(key,val)
                    self[key] = [val,unit]
                    self.write( "%s/%s.in"%(self.folder,filename) )
                    run(filename+".in")
                continue
            if type(values[0])==int:
                for val in values[1:]:
                    filename = "%s_%05d"%(key,val)
                    self[key] = [val,unit]
                    self.write( "%s/%s.in"%(self.folder,filename) )
                    run(filename+".in")
                continue
            if type(values[0])==list:
                for array in values[1:]:
                    filename = "%s_%s"%(key,"_".join(map(str,array)))
                    self[key] = [array,unit]
                    self.write( "%s/%s.in"%(self.folder,filename) )
                    run(filename+".in")
                continue
            if type(values[0])==complex:
                for value in values[1:]:
                    filename = "%s_%lf_%lf"%(key,value.real,value.imag)
                    self[key] = [value,unit]
                    self.write( "%s/%s.in"%(self.folder,filename) )
                    run(filename+".in")
                continue
            raise ValueError( "unknown type for variable:", key )

        #put back the original values of the variables
        for var in variables:
            self[var] = backup[var]

        return name_files

    def write(self,filename='yambo.in'):
        """
        Write a yambo input file
        """
        f = open(filename,"w")
        f.write(str(self))
        f.close()

    def pack(self,filename):
        """
        Pack all the data of this structure in a `.json` file
        """
        f = open(filename,'w')
        json.dump(f,[self.arguments,self.real,self.string,self.complex,self.array],indent=5)
        f.close()

    def __str__(self):
        """
        Returns the input file as a string
        """
        s  = ""

        #arguments
        s += "\n".join(self.arguments)+'\n'

        for key,value in list(self.variables.items()):
            if isinstance(value,basestring):
                s+= "%s = %10s\n"%(key,"'%s'"%value)
                continue
            if isinstance(value[0],float):
                val, unit = value
                s+="%s = %lf %s\n"%(key,val,unit)
                continue
            if isinstance(value[0],int):
                val, unit = value
                s+="%s = %d %s\n"%(key,val,unit)
                continue
            if isinstance(value[0],list):
                array, unit = value
                if isinstance(array[0],list):
                    s+='%% %s\n'%key
                    for l in array:
                        s+="%s \n"%(" | ".join(map(str,l))+' | ')
                    s+='%s'%unit
                    s+='%\n'
                else:
                    s+="%% %s\n %s %s \n%%\n"%(key," | ".join(map(str,array))+' | ',unit)
                continue
            if isinstance(value[0],basestring):
                array = value
                s+="%% %s\n %s \n%%\n"%(key," | ".join(["'%s'"%x.replace("'","").replace("\"","") for x in array])+' | ')
                continue
            if isinstance(value[0],complex):
                c, unit = value
                s+="%s = (%lf,%lf) %s\n"%(key,c.real,c.imag,unit)
                continue
            raise ValueError( "Unknown type %s for variable: %s" %( type(value), key) )
        return s
