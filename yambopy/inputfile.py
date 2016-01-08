# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
from subprocess import Popen, PIPE
import os
import json
from time import sleep
import re

class YamboIn():
    #Regular expressions
    _variaexp   = '([A-Za-z\_0-9]+(?:\_[A-Za-z]+)?)' #variables names 
    _numexp     = '([+-]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)' #number
    _spacexp    = '(?:\s+)?' #space
    _stringexp  = '"(.+)"' #string
    _arrayexp   = '%(?:\s+)?'+_variaexp+'\s+(?:\#.+)?((?:(?:\s|\.|[+-]?\d)+?\|)+)\s+([a-zA-Z]+)?' #arrays
    _complexexp = '\(\s+?'+_numexp+'\s+?,\s+?'+_numexp+'\s+?\)\s+([a-zA-Z]+)?' #complex numbers
    _runexp     = '([a-zA-Z0-9_]+)\s+' #runlevels
    # list of available runlevels to be stored in the arguments array
    _runlevels  = ['rim_cut','em1s','bse','optics','bsk','bss',
                   'em1d','gw0','HF_and_locXC','setup','ppa','cohsex','life',
                   'collisions','negf','el_ph_scatt','el_el_scatt','excitons','wavefunction','fixsymms'] 

    def __init__(self,args='',folder='.',vim=True,filename='yambo.in'):
        self.folder = folder

        #the type of the variables is determined from the type of variable in this dictionary
        self.variables = {} #here we will store the values of the variables
        self.arguments = [] #here we will store the arguments

        # if we initalize the class with arguments we call yambo to generate the input file
        if args != '':
            # if yambo calls vim we have to close it. We just want the generic input file
            # that yambo generates. Yambo should have a command line argument that just generates
            # the input file without calling the editor
            if vim:
                workdir = os.getcwd()
                os.chdir(folder)
                os.system('rm -f %s'%filename)
                yambo = Popen(args, stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True)
                yambo.stdin.write(":wq\n")
                yambo.stdin.flush()
                n = 0
                while not os.path.isfile(filename) and n < 6:
                    sleep(1.0)
                    n+=1
                yambo.kill()
                os.chdir(workdir)
            # if yambo is not compiled with vim we don't care
            else:
                os.system(args)
            self.read_file(filename="%s/yambo.in"%folder)

    def __getitem__(self,key):
        """ Get the value of a keyword in the input file
        """
        return self.variables[key]

    def __setitem__(self,key,value):
        """ Set the value of a keyword in the input file
        """
        self.variables[key] = value

    def __delitem__(self,key):
        """ remove a keyword from teh dicitonary 
        """
        del self.variables[key]    

    def read_file(self,filename='yambo.in'):
        """ Read the keywords from a file
        """
        yambofile = open(filename,"r")
        inputfile = self.read_string(yambofile.read())
        yambofile.close()

    def read_dict_type(self,variables):
        """ read the variables from a dictionary separated by datatype.
        This is done because the datatypes are not kept in 
        """
        self.variables = variables

    def add_dict(self,variables):
        """ add a dictionary to the current inputfile
        """
        self.variables.update(variables)

    def read_string(self,inputfile):
        """ Read the input variables from an input file
        """
        var_real     = re.findall(self._variaexp +'='+ self._spacexp +
                                  self._numexp + self._spacexp + '([A-Za-z]+)?',inputfile)  
        var_string   = re.findall(self._variaexp +'='+ self._spacexp + self._stringexp, inputfile)
        var_array    = re.findall(self._arrayexp, inputfile) 
        var_complex  = re.findall(self._variaexp +'='+ self._spacexp + self._complexexp, inputfile)
        var_runlevel = re.findall(self._runexp, inputfile)

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
            array = [val.strip() for val in array.split('|')[:-1]]
            self[name] = [array,unit]

        return {"arguments": self.arguments, "variables": self.variables}

    def optimize(self,conv,variables=('all',),run=lambda x: None):
        """ Function to to make multiple runs of yambo to converge calculation parameters
            Input:
            A dictionary conv that has all the variables to be optimized
            A list fo the name of the variables in the dicitonary that are to be optimized 
            A function run that takes as input the name of the inputfile (used to run yambo)

            Example:
            def run(filename):
                os.system('yambo -F %s'%filename)

            What is done
            Make a calculation with all the first elements of the different variables in conv
            For each of the remaining elements make a new run
        """
        name_files = []

        #check which variables to optimize
        if 'all' in variables:
            variables = conv.keys()

        #save all the variables
        backup = {}
        for var in variables:
            backup[var] = self[var]

        #add units to all the variables (to be improved!)
        for key,value in conv.items():
            if type(value[-1]) != str and type(value[0]) == list:
                conv[key] = [value,'']

        #make a first run with all the first elements
        reference = {}
        for key,value in conv.items():
            values, unit = value
            reference[key] = [values[0],unit]
            self[key] = [values[0],unit]
        #write the file and run
        self.write( "%s/reference.in"%(self.folder) )
        run('reference.in')

        #converge one by one
        for key in [var for var in conv.keys() if var in variables]:
            values, unit = conv[key]
            #put back the original values of the variables
            for var in variables:
                self[var] = reference[var]
            #change the other variables
            if type(values[0])==str:
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
            print "unknown type for variable:", key
            exit(1)

        #put back the original values of the variables
        for var in variables:
            self[var] = backup[var]

        return name_files 

    def write(self,filename='yambo.in'):
        """ write the a yambo input file 
        """
        f = open(filename,"w")
        f.write(str(self))
        f.close()    

    def pack(self,filename):
        """ pack all the data of this structure in a json file
        """
        f = open(filename,'w')
        json.dump(f,[self.arguments,self.real,self.string,self.complex,self.array],indent=5)
        f.close()

    def __str__(self):
        """  Function that returns the input file as a string
        """
        s  = ""

        #arguments
        s += "\n".join(self.arguments)+'\n'

        for key,value in self.variables.items():
            if type(value)==str:
                s+= "%s = %10s\n"%(key,"'%s'"%value)
                continue
            if type(value[0])==float:
                val, unit = value
                s+="%s = %lf %s\n"%(key,val,unit)
                continue
            if type(value[0])==int:
                val, unit = value
                s+="%s = %d %s\n"%(key,val,unit)
                continue
            if type(value[0])==list:
                array, unit = value
                s+="%% %s\n %s %s \n%%\n"%(key,"|".join(map(str,array))+'|',unit)
                continue
            if type(value[0])==complex:
                value, unit = value
                s+="%s = (%lf,%lf) %s\n"%(key,value.real,value.imag,unit)
                continue
            print "Unknown type for variable:", key
            exit(1)
        return s 
