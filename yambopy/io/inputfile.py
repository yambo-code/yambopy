# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
import os
import json
import re
from subprocess import Popen, PIPE
from yambopy import yambopyenv
from yambopy.tools.duck import isstring

def issave(path):
    """ Check if yambo SAVE folder is present either as directory or as symlink """
    if os.path.isdir(path) or os.path.islink(path): return True
    else: return False

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
    _runlevels  = ['rim_cut','chi','em1s','bse','optics','bsk','bss','dipoles','ExcitonGkkp'
                   'em1d','gw0','HF_and_locXC','setup','ppa','cohsex','life',
                   'collisions','negf','el_ph_scatt','el_el_scatt','excitons','wavefunction','fixsyms',
                   'QPDBs', 'QPDB_merge','RealTime','RT_X','RToccDos','RToccBnd','RToccEner',
                   'RToccTime','RTlifeBnd','amplitude','bzgrids','Random_Grid','gkkp','el_ph_corr','WRbsWF',
                   'Select_energy','RTDBs','photolum','kpts_map','RTtime','RToccupations','RTfitbands']

    def __init__(self,args=''):
        """
        Initalize the class
        """
        self.yamboargs = args

        #the type of the variables is determined from the type of variable in this dictionary
        self.variables = {} #here we will store the values of the variables
        self.arguments = [] #here we will store the arguments

    @classmethod
    def from_file(cls,filename='yambo.in',folder='.'):
        """ Create an instance of YamboIn from an existing input file"""
        instance = cls()
        err = instance.read_file(os.path.join(folder,filename))
        if err: raise FileNotFoundError('Could not read the %s file'%os.path.join(folder,filename))
        return instance

    @classmethod
    def from_runlevel(cls,runlevel,executable=yambopyenv.YAMBO,folder='.',filename='yambo.in'):
        """
        Create an input file from the runlevel.
        Will execute yambo in the folder with the runlevel arguments,
        read the file and return an instance of this class
        """
        workdir = os.getcwd()
        if executable[3:]=='ypp': filename='ypp.in'

        #check if there exists a SAVE folder
        save_path = os.path.join(folder,'SAVE')
        if not issave(save_path): raise ValueError('SAVE folder not found in %s'%save_path)

        #run yambo
        os.chdir(folder)
        if os.path.isfile(filename): os.remove(filename)
        if '-Q' not in runlevel: runlevel += ' -Q'
        if filename is not ('yambo.in' or 'ypp.in'): runlevel += ' -F %s'%filename
        command = "%s %s"%(executable,runlevel)
        yambo = Popen(command, stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True)
        yambo.wait()
        os.chdir(workdir)

        instance = cls()
        err = instance.read_file(os.path.join(folder,filename))

        if err:
            lines = []; app = lines.append
            app('Yambo did not create the %s input file.'%filename)
            app('command: %s'%command)
            app('folder:  %s/'%folder)
            raise FileNotFoundError("\n".join(lines))

        #read input file
        return cls.from_file(filename=filename,folder=folder)

    @classmethod
    def from_dictionary(cls,yamboin_dict):
        """Return an instance of this class from a dictionary"""
        yamboin = cls()
        yamboin.set_fromdict(yamboin_dict)
        return yamboin

    def set_fromdict(self,yamboin_dict):
        """Write a python script to generate this input"""
        #monkey patch the input file
        if 'arguments' in yamboin_dict.keys() and 'variables' in yamboin_dict.keys():
            for var,value in yamboin_dict.items():
                setattr(self,var,value)
        else:
            for var,value in yamboin_dict.items():
                self[var] = value

    def set_fromargs(self,yamboin_args):
        """Write a python script to generate this input"""
        #monkey patch the input file
        self.arguments += yamboin_args

    def __getitem__(self,key):
        """ Get the value of a variable in the input file
        """
        return self.variables[key]

    def __setitem__(self,key,value):
        """ Set the value of a variable in the input file
        """
        if key == "QPbands":
            #read from current QPkrange
            kstart,kstop,bstart,bstop = self.variables["QPkrange"][0]
            #read bands from QPbands
            new_bstart,new_bstop = value
            #set new QPkrange
            key = "QPkrange"
            value = [kstart,kstop,new_bstart,new_bstop]
        #if the units are not specified, add them
        if isinstance(value,list):
            if not any( [isinstance(v,str) for v in value] ):
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
        if not os.path.isfile(filename): return 1

        with open(filename,"r") as yambofile:
            inputfile = self.read_string(yambofile.read())
        return 0

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

    def optimize(self,conv,folder='.',variables=('all',),run=lambda x: None,ref_run=True):
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
            if not isinstance(value[-1],str) and type(value[0]) == list:
                conv[key] = [value,'']

        #make a first run with all the first elements
        reference = {}
        for key,value in list(conv.items()):
            values, unit = value
            reference[key] = [values[0],unit]
            self[key] = [values[0],unit]
        #write the file and run
        if ref_run==True:
            self.write( "%s/reference.in"%(folder) )
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
            if isinstance(values[0],str):
                for string in values[1:]:
                    filename = "%s_%s"%(key,string)
                    self[key] = string
                    self.write( folder + filename )
                    run(filename+".in")
                continue
            if type(values[0])==float:
                for val in values[1:]:
                    filename = "%s_%12.8lf"%(key,val)
                    self[key] = [val,unit]
                    self.write( "%s/%s.in"%(folder,filename) )
                    run(filename+".in")
                continue
            if type(values[0])==int:
                for val in values[1:]:
                    filename = "%s_%05d"%(key,val)
                    self[key] = [val,unit]
                    self.write( "%s/%s.in"%(folder,filename) )
                    run(filename+".in")
                continue
            if type(values[0])==list:
                for array in values[1:]:
                    filename = "%s_%s"%(key,"_".join(map(str,array)))
                    self[key] = [array,unit]
                    self.write( "%s/%s.in"%(folder,filename) )
                    run(filename+".in")
                continue
            if type(values[0])==complex:
                for value in values[1:]:
                    filename = "%s_%lf_%lf"%(key,value.real,value.imag)
                    self[key] = [value,unit]
                    self.write( "%s/%s.in"%(folder,filename) )
                    run(filename+".in")
                continue
            raise ValueError( "unknown type for variable:", key )

        #put back the original values of the variables
        for var in variables:
            self[var] = backup[var]

        return name_files

    def copy(self):
        """Return a copy of this object"""
        import copy
        return copy.deepcopy(self)

    def write(self,filename='yambo.in', prefix=''):
        """
        Write a yambo input file
        """
        with open(filename,"w") as f:
           f.write(prefix)
           f.write(str(self))

    def set_q(self,q):
        """Change one of ['QpntsRXp','QpntsRXd','QpntsRXs'] variables to calculate only one q-point"""
        for var in ['QpntsRXp','QpntsRXd','QpntsRXs']:
            qpts = self.variables.get(var,None)
            if qpts is None: continue
            self.variables[var] = [[q,q],'']
            return 0
        raise ValueError('Could not find one of the following variables set in the input file: \'QpntsRXp\',\'QpntsRXd\',\'QpntsRXs\'')

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
            if key == 'DrudeWXd':
                s+= 'DrudeWXd= '+value+"\n"
                continue
            if isstring(value):
                s+= "%s = %10s\n"%(key,"'%s'"%value)
                continue
            if isinstance(value[0],float):
                val, unit = value
                if val > 1e-6:
                    s+="%s = %lf %s\n"%(key,val,unit)
                else:
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
            if isinstance(value[0],str):
                array = value
                s+="%% %s\n %s \n%%\n"%(key," | ".join(["'%s'"%x.replace("'","").replace("\"","") for x in array])+' | ')
                continue
            if isinstance(value[0],complex):
                c, unit = value
                s+="%s = (%lf,%lf) %s\n"%(key,c.real,c.imag,unit)
                continue
            raise ValueError( "Unknown type %s for variable: %s" %( type(value), key) )
        return s
