"""
Module that manages the parsing of the ``ndb.RT_carriers`` database created by `yambo_rt`.
"""
from netCDF4 import Dataset
from yambopy.units import ha2ev
from yambopy.plot.spectra import get_spectra

import numpy as np

class YamboRT_Carriers_DB():
    """
    Class to manage information about the real time distrubtion of carriers from the
    ``ndb.RT_carriers`` database created by `yambo_rt`.

    Args:
        file (:py:class:`string`) : string with the name of the database to be parsed
        verbose (:py:class:`boolean`) : define the amount of information provided on terminal

    Attributes:
        E_bare (:py:class:`np.array`) : Array that contains the bare bands energies.
            The structure of the array is [numk*numbands]. Energies are expressed in eV
            The structure is k1*b1,k1*b2,..,k1*bn,k2*b1,...
        f_bare (:py:class:`np.array`) : Array that contains the bare bands occupations.
            The structure of the array is [numk*numbands]
        kpoints (:py:class:`np.array`) : Array that contains the kpoints in cartesian  coordinates
            in units of 2*np.pi/alat (with a vector alat). The structure of the array is [numk,3],
            the first index runs over the kpoints and the second one gives the component
        bands_kpts (:py:class:`np.array`) : Array with the [first_band,last_band,numk]
        k_weight (:py:class:`np.array`) : Array that contains the weights of the kpoints.
            The structure of the array is [numk]
        delta_E (:py:class:`np.array`) : Array that contains the time-dependent correction
            to the bands energies. The structure of the array is [time,numk*numbands].
            Energies are expressed in eV
        delta_f (:py:class:`np.array`) : Array that contains the time-dependent corrections to
            the bands occupations. The structure of the array is [time,numk*numbands]

    """

    def __init__(self,folder='.',calc='SAVE',carriers_db='ndb.RT_carriers'):
        # Find path with RT data
        self.carriers_path = '%s/%s/%s'%(folder,calc,carriers_db)
        self.calc=calc
        try:
            data_obs= Dataset(self.carriers_path)
        except:
            raise ValueError("Error reading CARRIERS database at %s"%self.carriers_path)

        self.readDB(data_obs)

        data_obs.close()


    def readDB(self,database):
        """
        Read the data from the ``ndb.RT_carriers`` database created by `yambo_rt`. The variables
        are extracted from the database and stored in the attributes of the object.

        Args:
            verbose (:py:class:`boolean`) : define the amount of information provided on terminal

        """

        self.E_bare = ha2ev*np.array(database.variables['RT_carriers_E_bare'])
        self.f_bare = np.array(database.variables['RT_carriers_f_bare'])
        self.kpoints = np.array(database.variables['RT_kpt'][:].T)
        self.bands_kpts = np.array(database.variables['RT_bands_kpts'])
        self.k_weight = np.array(database.variables['RT_k_weight'])
        self.delta_E = ha2ev*np.array(database.variables['RT_carriers_delta_E'])
        self.delta_f = np.hstack(np.array(database.variables['RT_carriers_delta_f']))

    def get_info(self):
        """
        Provide information on the attributes of the class
        """
        print('YamboRTCarriersParser variables structure')
        print('Bands used and number of k-points',self.bands_kpts)
        print('kpoints shape',self.kpoints.shape)
        print('E_bare shape',self.E_bare.shape)
        print('f_bare shape',self.f_bare.shape)
        print('delta_E shape',self.delta_E.shape)
        print('delta_f',self.delta_f.shape)


    def build_f_bare_dos(self, dE = 0.1, eta = 0.1, broad_kind = 'lorentzian'):
        """
        For each kpoint build a dos which expresses the bare occupation level in terms
        of the energy. The energy ranges from the minum to the maximum of the
        E_bare variable.

        Args:
            dE (:py:class:`float`) : energy step in eV
            eta (:py:class:`float`) : magnitude of the broading parameter (in the same units used for the values array)
            broad_kind (:py:class:`string`) : type of broading function used (lorentzian, gaussian)

        Returns:
            :py:class:`Dos` : Instance of the ``Dos`` class. The object is an array of dos, one for
                each kpoint

        """
        numbnds = self.bands_kpts[1]-self.bands_kpts[0]
        numkp = self.bands_kpts[2]
        Emin= min(self.E_bare)-10*eta
        Emax= max(max(self.E_bare*self.f_bare),max(self.E_bare*self.delta_f))+10*eta
        #
        # Here I want to calculate the weighted density of states (w-DOS) so I pass
        # to the get_spectra function the occupation as residual for the DOS
        #
        w, dos = get_spectra(energies=self.E_bare,residuals=self.f_bare,broadening=eta,emin=Emin,emax=Emax)
        return w,dos
    
    
    def build_delta_f_dos(self, dE = 0.1, eta = 0.1, broad_kind = 'lorentzian'):
        """
        For each kpoint build a dos which expresses the bare occupation level in terms
        of the energy. The energy ranges from the minum to the maximum of the
        E_bare variable.

        Args:
            dE (:py:class:`float`) : energy step in eV
            eta (:py:class:`float`) : magnitude of the broading parameter (in the same units used for the values array)
            broad_kind (:py:class:`string`) : type of broading function used (lorentzian, gaussian)

        Returns:
            :py:class:`Dos` : Instance of the ``Dos`` class. The object is an array of dos, one for
                each kpoint

        """
        numbnds = self.bands_kpts[1]-self.bands_kpts[0]
        numkp = self.bands_kpts[2]
        Emin= min(self.E_bare)-10*eta
        Emax= max(max(self.E_bare*self.f_bare),max(self.E_bare*self.delta_f))+10*eta
        #
        # Here I want to calculate the weighted density of states (w-DOS) so I pass
        # to the get_spectra function the occupation as residual for the DOS
        #
        w, dos = get_spectra(energies=self.E_bare,residuals=self.delta_f,broadening=eta,emin=Emin,emax=Emax)
        return w,dos


