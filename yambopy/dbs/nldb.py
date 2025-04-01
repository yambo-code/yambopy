# Copyright (c) 2023, Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.plot import *
from yambopy.units import ha2ev,fs2aut,speed_of_light
import numpy as np
import sys
import os

#
# This class reads all data from the ndb.Nonlinear database
# and its fragment ndb.Nonlinear_fragment_xxx
# 
# Time series is stored in IO_TIME_points variable
# while external fields, current and polarization of the different runs 
# are stored in Efield[x],Polarization[x],Current[x]
# It is a single run just use Efield[0],Current[0] etc...
#
class YamboNLDB(object):
    """
    Open the NL databases and store it in a NLDB class.
    """

    def __init__(self,folder='.',calc='SAVE',nl_db='ndb.Nonlinear'):
        # Find path with RT data
        self.nl_path = '%s/%s/%s'%(folder,calc,nl_db)
        self.calc=calc
        try:
            data_obs= Dataset(self.nl_path)
        except:
            raise ValueError("Error reading OBSERVABLES database at %s"%self.nl_path)

        self.read_observables(data_obs)


        data_obs.close()

    def read_Efield(self,database,RT_step,n):
         efield={}
         efield["name"]       =database.variables['Field_Name_'+str(n)][...].tostring().decode().strip()
         efield["versor"]     =database.variables['Field_Versor_'+str(n)][:].astype(np.double)
         efield["intensity"]  =database.variables['Field_Intensity_'+str(n)][0].astype(np.double)
         efield["damping"]    =database.variables['Field_Damping_'+str(n)][0].astype(np.double)
         efield["freq_range"] =database.variables['Field_Freq_range_'+str(n)][:].astype(np.double)
         efield["freq_steps"] =database.variables['Field_Freq_steps_'+str(n)][:].astype(np.double)
         efield["freq_step"]  =database.variables['Field_Freq_step_'+str(n)][0].astype(np.double)
         efield["initial_time"]  =database.variables['Field_Initial_time_'+str(n)][0].astype(np.double)
         #
         # set t_initial according to Yambo 
         #
         efield["initial_indx"] =max(round(efield["initial_time"]/RT_step)+1,2)
         efield["initial_time"] =(efield["initial_indx"]-1)*RT_step
         #
         # define the field amplitude
         #
         efield["amplitude"]    =np.sqrt(efield["intensity"]*4.0*np.pi/speed_of_light)

         return efield

    def read_observables(self,database):
        """
        Read all data from the database
        """
        self.Gauge          = database.variables['GAUGE'][...].tostring().decode().strip()
        self.NE_steps       = database.variables['NE_steps'][0].astype('int')
        self.RT_step        = database.variables['RT_step'][0].astype(np.double)
        self.n_frequencies  = database.variables['n_frequencies'][0].astype('int')

        try:
            self.n_angles       = database.variables['n_angles'][0].astype('int')
        except:
            self.n_angles   = 0
        try:
            self.NL_initial_versor = database.variables['NL_initial_versor'][:].astype(np.double)
        except:
            self.NL_initial_versor = [0.0, 0.0, 0.0] 

        self.NL_damping     = database.variables['NL_damping'][0].astype(np.double)
        self.RT_bands       = database.variables['RT_bands'][:].astype('int')
        self.NL_er          = database.variables['NL_er'][:].astype(np.double)
        self.l_force_SndOrd = database.variables['l_force_SndOrd'][0].astype('bool')
        self.l_use_DIPOLES  = database.variables['l_use_DIPOLES'][0].astype('bool')
        try:
            self.l_eval_CURRENT = database.variables['l_eval_CURRENT'][0].astype('bool')
        except:
            self.l_eval_CURRENT = False
        self.QP_ng_SH       = database.variables['QP_ng_SH'][0].astype('int')
        self.QP_ng_Sx       = database.variables['QP_ng_Sx'][0].astype('int')
        self.RAD_LifeTime   = database.variables['RAD_LifeTime'][0].astype(np.double)
        self.Integrator     = database.variables['Integrator'][...].tostring().decode().strip()
        self.Correlation    = database.variables['Correlation'][...].tostring().decode().strip()
        #
        # Time variables
        #
        self.IO_TIME_N_points  = database.variables['IO_TIME_N_points'][0].astype('int')
        self.IO_TIME_LAST_POINT= database.variables['IO_TIME_LAST_POINT'][0].astype('int')
        self.IO_TIME_points    = database.variables['IO_TIME_points'][:].astype(np.double)
        #
        # External fields
        # 
        self.Efield_general=[]
        for n in range(1,4):
            efield=self.read_Efield(database,self.RT_step,n)
            self.Efield_general.append(efield.copy())

        #
        # Read polarization and currect files 
        #
        self.Polarization=[]
        self.Current     =[]
        self.E_ext       =[]
        self.E_tot       =[]
        self.E_ks        =[]
        self.Efield      =[] # Store the first external field for each run at different frequencies
        self.Efield2     =[]
        #
        if self.n_angles!=0:
            self.n_runs=self.n_angles
        if self.n_frequencies!=0:
            self.n_runs=self.n_frequencies
        if (self.n_angles!=0 and self.n_frequencies!=0):
            print("Error both n_angles and n_frequencies !=0 ")
            sys.exit(0)
        if (self.n_angles==0 and self.n_frequencies==0):
            self.n_runs=1

        #
        for f in range(self.n_runs):
            try:
                data_p_and_j= Dataset(self.nl_path+"_fragment_"+str(f+1))
            except:
                print("Error reading database: %s" % self.nl_path+"_fragment_"+str(f+1))
                continue
            pol  = data_p_and_j.variables['NL_P_freq_'+str(f+1).zfill(4)][:,:].astype(np.double)
            curr = data_p_and_j.variables['NL_J_freq_'+str(f+1).zfill(4)][:,:].astype(np.double)
            e_ext= data_p_and_j.variables['E_ext_freq_'+str(f+1).zfill(4)][:,:,:].astype(np.double)
            e_ext_c=e_ext[:,:,0]+1j*e_ext[:,:,1]
            e_tot= data_p_and_j.variables['E_tot_freq_'+str(f+1).zfill(4)][:,:,:].astype(np.double)
            e_tot_c=e_tot[:,:,0]+1j*e_tot[:,:,1]
            e_ks = data_p_and_j.variables['E_ks_freq_'+str(f+1).zfill(4)][:,:,:].astype(np.double)
            e_ks_c=e_ks[:,:,0]+1j*e_ks[:,:,1]

            self.Polarization.append(pol.copy())
            self.Current.append(curr.copy())
            self.E_ext.append(e_ext_c.copy())
            self.E_tot.append(e_tot_c.copy())
            self.E_ks.append(e_ks_c.copy())

            # Read only the first field for SHG
            # I don't need it in the pump-probe configuration
            efield=self.read_Efield(data_p_and_j,self.RT_step,1)
            efield2=self.read_Efield(data_p_and_j,self.RT_step,2)
            self.Efield.append(efield.copy())
            self.Efield2.append(efield2.copy())

    def __str__(self):
        """
        Print all info of the database
        """
        s="\n * * * ndb.Nonlinear db data * * * \n\n"
        s+="Gauge         : "+str(self.Gauge)+"\n"
        s+="NE_steps      : "+str(self.NE_steps)+"\n"
        s+="RT_step       : "+str(self.RT_step/fs2aut)+" [fs] \n"
        s+="n_frequencies : "+str(self.n_frequencies)+"\n"   
        s+="n_angles      : "+str(self.n_angles)+"\n"   
        s+="NL_initial_versor   : "+str(self.NL_initial_versor)+"\n"   
        s+="NL_damping    : "+str(self.NL_damping*ha2ev)+"\n"
        s+="RT_bands      : "+str(self.RT_bands)+"\n"
        s+="NL_er         : "+str(self.NL_er*ha2ev)+" [eV] \n"
        s+="Sencond Order : "+str(self.l_force_SndOrd)+"\n" 
        s+="Use Dipoles   : "+str(self.l_use_DIPOLES)+"\n"
        s+="QP_ng_SH      : "+str(self.QP_ng_SH)+"\n"
        s+="QP_ng_Sx      : "+str(self.QP_ng_Sx)+"\n"  
        s+="RAD_LifeTime  : "+str(self.RAD_LifeTime/fs2aut)+" [fs] \n" 
        s+="Integrator    : "+str(self.Integrator)+"\n"
        s+="Correlation   : "+str(self.Correlation)+"\n"
        for efield in self.Efield:
            if efield["name"] == "none":
                continue
            s+="\nEfield name         : "+str(efield["name"])+"\n"
            s+="Efield versor       : "+str(efield["versor"])+"\n"  
            s+="Efield Intesity     : "+str(efield["intensity"])+"\n"
            s+="Efield Damping      : "+str(efield["damping"])+"\n"
            s+="Efield Freq range   : "+str(efield["freq_range"]*ha2ev)+" [ev] \n"
            s+="Efield Initial time : "+str(efield["initial_time"]/fs2aut)+" [fs] \n"
        return s

