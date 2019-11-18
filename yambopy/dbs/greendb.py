# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
import shutil
ha2ev  = 27.211396132

class YamboGreenDB(object):
    """
    Read the green's functions calculated using yambo
    These green's functions describe the spectral function of the quasiparticles.
    The quasi-particles can be from electron-phonon or GW calculations
    """
    def __init__(self,filename='ndb.G',folder='.'):
        self.folder = folder
        self.filename = "%s/%s"%(folder,filename)

        #read em1s database
        try:
            database = Dataset(self.filename, 'r')
        except:
            raise IOError("Error opening %s in YamboGreenDB"%self.filename)

        #read the Green's functions energies
        re,im = database.variables['Green_Functions_Energies'][:]
        self.energies = (re+im*1j).T

        #read the Green's Functions
        re,im = database.variables['Green_Functions'][:]
        self.green = (re+im*1j).T
        
        #read the self-energy operator
        re,im = database.variables['SE_Operator'][:]
        self.se = (re+im*1j).T

        self.nqps, self.nenergies = self.green.shape

        #read QP_table
        qptable = database.variables['QP_table'][:].astype(int)
        self.band1, self.band2, self.kindex = qptable
        self.bandmax = max(self.band1)
        self.bandmin = min(self.band1)

        #qp dictionary
        self.qp_dict = {}
        for nqp,(b1,b2,kindex) in enumerate(qptable.T):
            self.qp_dict[(b1,b2,kindex)] = nqp

        #read QP_kpts
        kpts = database.variables['QP_kpts'][:].T
        self.qpoints = kpts.shape

    def plot(self,ax,kpt=0,band=0,what='SE',e0=None,**kwargs):
        """
        Plot quantities from this database
        """
        nqp = self.qp_dict[(band,band,kpt)]
        
        x = self.energies[nqp]
        options = {'SE':self.se,
                   'green':self.green}
        y = options[what][nqp]

        #get band and k-point
        band = self.band1[nqp]
        kpt  = self.kindex[nqp]

        ax.set_title('kpt=%d band=%d'%(kpt,band))
        ax.plot(x.real,y.real,label='Re(%s)'%what,**kwargs)
        ax.plot(x.real,y.imag,label='Im(%s)'%what,**kwargs)
        if e0 is not None:
            ax.plot(x.real,e0[nqp]-x.real)

            #plot 0
            ax.axhline(0,c='k',lw=1)

            #set axis
            rmin, rmax = min(y.real),max(y.real)
            imin, imax = min(y.imag),max(y.imag)
            ax.set_ylim(min(rmin,imin),max(rmax,imax))

    def modQP(self,filename_reference,filename_new):
        """
        Take a QP file as reference and modify the values of the energies, lifetimes and Z factors
        according to the ones calculated from ndb.Green.
        
        Arguments:
            filename_reference : name of the reference file
            filename_new : name of the new file with the calculated data
        """

        #copy ref file to new file
        shutil.copy(filename_reference, filename_new)

        #read QP file
        qp = Dataset(filename_new,'r+')

        #check dimensions
        #print qp.variables['QP_E_Eo_Z'][:].shape
        #print self.eqp.shape
        #print self.z.shape

        qp.variables['QP_E_Eo_Z'][0,:,0] = self.eqp.real
        qp.variables['QP_E_Eo_Z'][1,:,0] = self.eqp.imag
        qp.variables['QP_E_Eo_Z'][0,:,2] = self.z.real
        qp.variables['QP_E_Eo_Z'][1,:,2] = self.z.imag

        #write 
        qp.close()

    def getQP(self,e0,bandmin=None,bandmax=None,debug=False,secant=True,braket=None):
        """
        Get quasiparticle states
    
        Arguments:
        e0 -> bare eigenvalues in eV
        """
        from scipy.optimize import bisect, newton
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
            x = np.linspace(emin,emax,100)
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

        if bandmin is None: bandmin = self.bandmin
        if bandmax is None: bandmax = self.bandmax
 
        self.eqp = np.zeros([self.nqps],dtype=complex) 
        self.z   = np.zeros([self.nqps],dtype=complex) 
        for nqp in range(self.nqps):

            band = self.band1[nqp]
            kpt  = self.kindex[nqp]
            if debug: print("%3d %3d %3d %8.4lf"%(nqp, kpt, band, e0[nqp]))

            if not (bandmin <= band <= bandmax):
                continue

            #get x and y
            x = self.energies[nqp].real
            y = self.se[nqp]

            #interpolate real part of function
            f = interp1d(x,y.real-x+e0[nqp],kind='slinear')

            #find zero
            if secant:
                try:
                    eqp = newton(f,e0[nqp],maxiter=200)
                except ValueError as msg:
                    print(msg)
                    if debug: error(nqp)
            else:
                if braket:
                    emin = e0[nqp]-braket
                    emax = e0[nqp]+braket
                else:
                    emin = min(x)
                    emax = max(x)

                eqp = bisect(f,emin,emax) 

            #interpolate whole function
            f = interp1d(x,y)

            #calculate Z factors
            dse = derivative(f,eqp,dx=1e-8)
            z = 1./(1-dse)

            #find Im(Se(EQP)) which corresponds to the lifetime
            lif = f(eqp).imag
            eqp += 1j*lif

            #store values
            self.eqp[nqp] = eqp
            self.z[nqp] = z           

            #cehck for potential errors
            if z>1 and debug:
                print(z)
                error(nqp)            

        return self.eqp, self.z

    def __str__(self):
        s = ""
        s += "nenergies: %d\n"%self.nenergies
        s += "nqps:      %d\n"%self.nqps
        s += "bandmin:   %d\n"%self.bandmin
        s += "bandmax:   %d"%self.bandmax
        return s

