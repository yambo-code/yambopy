from __future__ import division
# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from builtins import range
from builtins import object
from past.utils import old_div
from yambopy import *

ha2ev  = 27.211396132

class YamboGreenDB(object):
    """
    Read the green's functions calculated using yambo
    These green's functions describe the spectral function of the quasiparticles.
    The quasi-particles can be from electron-phonon or GW calculations
    """
    def __init__(self,save='SAVE',filename='ndb.G'):
        self.filename = "%s/%s"%(save,filename)

        #read em1s database
        try:
            database = Dataset(self.filename, 'r')
        except:
            raise IOError("Error opening %s in YamboGreenDB"%self.filename)

        #read the Green's functions energies
        re,im = database['Green_Functions_Energies'][:]*ha2ev
        self.energies = (re+im*1j).T

        #read the Green's Functions
        re,im = database['Green_Functions'][:]*ha2ev
        self.green = (re+im*1j).T
        
        #read the self-energy operator
        re,im = database['SE_Operator'][:]*ha2ev
        self.se = (re+im*1j).T

        self.nqps, self.nenergies = self.green.shape

        #read QP_table
        qptable = database['QP_table'][:]
        self.band1, self.band2, self.kindex = qptable

        #read QP_kpts
        kpts = database['QP_kpts'][:].T
        self.qpoints = kpts.shape

    def plot(self,ax,nqp=0,nb=0,what='SE',**kwargs):
        """
        Plot quantities from this database
        """
        x = self.energies[nqp]
        options = {'SE':self.se,
                   'green':self.green}
        y = options[what][nqp]

        ax.plot(x.real,y.real,label='Re(%s)'%what,**kwargs)
        ax.plot(x.real,y.imag,label='Im(%s)'%what,**kwargs)

    def getQP(self,e0,debug=False):
        """
        Get quasiparticle states
    
        Arguments:
        e0 -> bare eigenvalues
        """
        from scipy.optimize import bisect
        from scipy.interpolate import interp1d
        from scipy.misc import derivative

        #check if the eigenvalues have the correct dimensions
        if len(e0) != self.nqps:
            raise ValueError('Wrong dimensions in bare eigenvalues')

        #in case something is strange we plot the stuff
        def error(nqp):
            ax = plt.gca()

            #plot 0
            ax.axhline(0,c='k',lw=1)
            
            #se limits
            semin = min(self.se[nqp].real)
            semax = max(self.se[nqp].real)
            plt.ylim(semin,semax)

            #plot self energy
            self.plot(ax,nqp=nqp)

            #plot omega-e0
            emin = min(self.energies[nqp].real)
            emax = max(self.energies[nqp].real)
            x = np.arange(emin,emax,old_div((emax-emin),100))
            plt.plot(x,x-e0[nqp])

            #plot imaginary part of greens funciton
            x = self.energies[nqp].real
            y = self.green[nqp].imag
            plt.plot(x,y/max(y)*semax)
    
            #plot eqp
            #plt.axvline(self.eqp[nqp],lw=1)
            #plt.axvline(e0[nqp],lw=1)

            plt.legend(frameon=False)
            plt.show()
 
        self.eqp = np.zeros([self.nqps],dtype=complex) 
        self.z   = np.zeros([self.nqps],dtype=complex) 
        for nqp in range(self.nqps):

            #get x and y
            x = self.energies[nqp].real
            y = self.se[nqp]

            #interpolate real part of function
            f = interp1d(x,y.real-x+e0[nqp],kind='slinear')

            #find zero
            eqp = bisect(f,min(x),max(x)) 

            #interpolate whole function
            f = interp1d(x,y)

            #calculate Z factors
            #Z = (1-dSE/de)^(-1)
            dse = derivative(f,eqp,dx=1e-8)
            z = old_div(1.,(1-dse))

            #find Im(Se(EQP)) which corresponds to the lifetime
            lif = f(eqp).imag
            eqp += 1j*lif

            #store values
            self.eqp[nqp] = eqp
            self.z[nqp] = z           

            #cehck for potential errors
            if z>1 and debug:
                error(nqp)            

        return self.eqp, self.z
 
    def __str__(self):
        s = ""
        s += "nenergies: %d\n"%self.nenergies
        s += "nqps:      %d"%self.nqps
        return s

