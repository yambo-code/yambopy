# Copyright (c) 2022. Yambopy team
# All rights reserved.
#
# This file is part of the yambopy project
# Authors: H. Miranda, A. Molina, F. Paleari
#
import os
import numpy as np
from yambopy.units import ha2ev
from yambopy.tools.string import marquee
from yambopy.plot.bandstructure import YambopyBandStructure, YambopyBandStructureList
from yambopy.plot.plotting import add_fig_kwargs
from yamboparser import YamboFile
from yambopy.lattice import car_red, red_car, rec_lat, vol_lat

class YamboQPDB():
    """
    Class to read yambo ndb.QP files

    These files describe the quasiparticle states calculated from yambo
    Includes the quasi-particle energies, the lifetimes and the Z factors
    """
    def __init__(self,qps):
        """
        Initialize the YamboQP class
        """
        self.qps          = qps
        self.kpoints_iku  = np.array(qps['Kpoint'])
        self.kpoint_index = np.array(qps['Kpoint_index'],dtype=int)
        self.band_index   = np.array(qps['Band'],dtype=int)
        self.e0           = np.array(qps['Eo']).real*ha2ev
        self.e            = np.array(qps['E']).real*ha2ev
        self.linewidths   = np.array(qps['E']).imag*ha2ev
        self.qpz          = np.array(qps['Z']).real
        self.spin         = False
        if 'Spin_pol' in qps.keys():
            self.spin_pol = qps['Spin_pol']
            self.spin     = True

    @property
    def eigenvalues_qp(self):
        if not hasattr(self,'_eigenvalues_qp'):
            self._eigenvalues_dft, self._eigenvalues_qp, self._lifetimes, self._z = self.get_qps()
        return self._eigenvalues_qp

    @property
    def eigenvalues_dft(self):
        if not hasattr(self,'_eigenvalues_dft'):
            self._eigenvalues_dft, self._eigenvalues_qp, self._lifetimes, self._z = self.get_qps()
        return self._eigenvalues_dft

    @property
    def lifetimes(self):
        if not hasattr(self,'_lifetimes'):
            self._eigenvalues_dft, self._eigenvalues_qp, self._lifetimes, self._z = self.get_qps()
        return self._lifetimes

    @property
    def z(self):
        if not hasattr(self,'_z'):
            self._eigenvalues_dft, self._eigenvalues_qp, self._lifetimes, self._z = self.get_qps()
        return self._z

    @classmethod
    def from_db(cls,filename='ndb.QP',folder='.'):
        """
        Create instance of this class from a ndb.QP file
        """
        cls.qp_path = folder
        db_path = os.path.join(folder,filename)
        if os.path.isfile(db_path):
            yfile = YamboFile(filename,folder)
        else:
            raise IOError('File %s not found'%db_path)
        return cls(yfile.data)
    from_db_file = from_db
    
    def get_qps(self):
        """
        Get quasiparticle energies in a list
        """
        #start arrays

        # AMS: I changed the way we define the arrays. Hope is not breaking other things

        # I have shifted 
        ncalculatedkpoints = self.max_kpoint - self.min_kpoint + 1
        if self.spin is True:
           eigenvalues_dft = np.zeros([ncalculatedkpoints,self.nbands,2])
           eigenvalues_qp  = np.zeros([ncalculatedkpoints,self.nbands,2])
           linewidths      = np.zeros([ncalculatedkpoints,self.nbands,2])
           z               = np.zeros([self.nkpoints,self.nbands,2])

           for ei,e0i,li,zi,ki,ni,s_pol in zip(self.e,self.e0,self.linewidths,self.qpz,self.kpoint_index,self.band_index,self.spin_pol):
                # position in array
                nkpoint = ki - self.min_kpoint 
                nband   = ni - self.min_band
                spol    = int(s_pol - 1)

                eigenvalues_dft[nkpoint,nband,spol] = e0i
                eigenvalues_qp[nkpoint,nband,spol]  = ei
                linewidths[nkpoint,nband,spol]      = li
                z[nkpoint,nband,spol]               = zi

        else:
            eigenvalues_dft = np.zeros([ncalculatedkpoints,self.nbands])
            eigenvalues_qp  = np.zeros([ncalculatedkpoints,self.nbands])
            linewidths      = np.zeros([ncalculatedkpoints,self.nbands])
            z               = np.zeros([self.nkpoints,self.nbands])

            for ei,e0i,li,zi,ki,ni in zip(self.e,self.e0,self.linewidths,self.qpz,self.kpoint_index,self.band_index):

                # position in array
                nkpoint = ki - self.min_kpoint 
                nband   = ni - self.min_band

                eigenvalues_dft[nkpoint,nband] = e0i
                eigenvalues_qp[nkpoint,nband]  = ei
                linewidths[nkpoint,nband]      = li
                z[nkpoint,nband]               = zi
        return eigenvalues_dft, eigenvalues_qp, linewidths, z


    def get_filtered_qps(self,min_band=None,max_band=None):
        """Return selected QP energies as a flat list"""
        e0=[]; qp=[]; lw=[]
        for ei,e0i,li,ki,ni in zip(self.e,self.e0,self.linewidths,self.kpoint_index,self.band_index):
            if min_band and ni < min_band: continue
            if max_band and ni > max_band: continue
            e0.append(e0i)
            qp.append(ei)
            lw.append(lw)
        return e0,qp,lw

    def get_direct_gaps(self,valence):
        """
        Compute the QP and DFT gaps
        
        Arguments:
            valence: number of bands in the valence
        """
        na = np.newaxis
        shifted_valence = valence-self.min_band+1
 
        #direct gap
        dft_jdos = self.eigenvalues_dft[:,na,shifted_valence:]-self.eigenvalues_dft[:,:shifted_valence,na]
        qp_jdos  = self.eigenvalues_qp[:,na,shifted_valence:] -self.eigenvalues_qp[:,:shifted_valence,na]
        direct_dft_gap = np.min(dft_jdos)
        direct_qp_gap  = np.min(qp_jdos)

        #indirect gap
        #TODO take the min and max of the VBM and CBM
        return direct_dft_gap, direct_qp_gap

    def get_scissor(self,valence,verbose=1):
        """
        Compute the scissor operator replicating the QP corrections
        
        Arguments:
            valence: number of bands in the valence
        """
        from scipy import stats
        lines = []; app = lines.append
        #valence
        e0, eqp, lw = self.get_filtered_qps(self.min_band,valence)
        vslope, vintercept, r_value, p_value, std_err = stats.linregress(e0,eqp)
        app('valence bands:')
        app('slope:     {}'.format(vslope))
        app('intercept: {}'.format(vintercept))
        app('r_value:   {}'.format(r_value))
        app('p_value:   {}'.format(p_value))
        app('std_err:   {}'.format(std_err))
       
        #conduction
        e0, eqp, lw = self.get_filtered_qps(valence+1,self.max_band)
        cslope, cintercept, r_value, p_value, std_err = stats.linregress(e0,eqp)
        app('\nconduction bands:')
        app('slope:     {}'.format(cslope))
        app('intercept: {}'.format(cintercept))
        app('r_value:   {}'.format(r_value))
        app('p_value:   {}'.format(p_value))
        app('std_err:   {}'.format(std_err))

        #get gaps
        direct_dft_gap,direct_qp_gap = self.get_direct_gaps(valence) 
        shift = direct_qp_gap-direct_dft_gap
        app('direct dft gap: {}'.format(direct_dft_gap))
        app('direct qp gap:  {}'.format(direct_qp_gap))

        scissor_list = [shift,cslope,vslope]
        app('\vscissor list (shift,c,v) [eV,adim,adim]: {}'.format(scissor_list))
        if verbose: print("\n".join(lines))
        return shift,cslope,vslope,cintercept,vintercept

    def plot_scissor_ax(self,ax,valence,verbose=1):
        # Create option to plot correction vs ks eigenvalues
        """
        Plot the scissor on a matplotlib axis
        """
        shift,cslope,vslope,cintercept,vintercept=self.get_scissor(valence,verbose=verbose)

        #plot qps
        ve0,vqp,_=self.get_filtered_qps(self.min_band,valence)
        ve0 = np.array(ve0)
        vqp = np.array(vqp)
        #ax.scatter(ve0,vqp)
        ce0,cqp,_=self.get_filtered_qps(valence+1,self.max_band)
        ce0 = np.array(ce0)
        cqp = np.array(cqp)
        #ax.scatter(ce0,cqp)

        #plot the fits
        vx = np.linspace(np.min(ve0),np.max(ve0),2)
        cx = np.linspace(np.min(ce0),np.max(ce0),2)
        vy = vslope*vx+vintercept
        cy = cslope*cx+cintercept
        #ax.plot(vx,vy)
        #ax.plot(cx,cy)

        return ax.scatter(ve0,vqp), ax.scatter(ce0,cqp), ax.plot(vx,vy), ax.plot(cx,cy)

    @add_fig_kwargs
    def plot_scissor(self,valence,verbose=1):
        """
        Plot the the QP energies and the scissor fit
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_scissor_ax(ax,valence,verbose=verbose)
        return fig

    def get_bs_path(self,lat,path,debug=False,**kwargs):
        """Get a band-structure on a path"""
        bands_kpoints, bands_indexes, path_car = lat.get_path(path.kpoints,debug=debug)

        # set fermi energy
        # NOT EVIDENT IN SPIN-POLARIZED SYSTEM
        #print(lat.electrons)
        #self.top_valence_band = int(lat.electrons/2)
        #print(self.top_valence_band)
        #print(self.top_valence_band-7)
        #print(self.eigenvalues_dft[:,0:(self.top_valence_band-7+2),0])
        #print(self.eigenvalues_dft[:,0:(self.top_valence_band-7+2),1])
        #exit()

        red_bands_kpoints = car_red(bands_kpoints,lat.rlat)
        if self.spin == True:
           print('Spin polarized bands')
           ks_bandstructure_up = YambopyBandStructure(self.eigenvalues_dft[bands_indexes,:,0],red_bands_kpoints,kpath=path,**kwargs)
           ks_bandstructure_dw = YambopyBandStructure(self.eigenvalues_dft[bands_indexes,:,1],red_bands_kpoints,kpath=path,**kwargs)
           qp_bandstructure_up = YambopyBandStructure(self.eigenvalues_qp[bands_indexes,:,0], red_bands_kpoints,kpath=path,**kwargs)
           qp_bandstructure_dw = YambopyBandStructure(self.eigenvalues_qp[bands_indexes,:,1], red_bands_kpoints,kpath=path,**kwargs)

           return ks_bandstructure_up, ks_bandstructure_dw, qp_bandstructure_up, qp_bandstructure_dw
        else:
           print('no polarized bands')
           ks_bandstructure = YambopyBandStructure(self.eigenvalues_dft[bands_indexes],red_bands_kpoints,kpath=path,**kwargs)
           qp_bandstructure = YambopyBandStructure(self.eigenvalues_qp[bands_indexes] ,red_bands_kpoints,kpath=path,**kwargs)

           return ks_bandstructure, qp_bandstructure
        
    def get_bs(self,**kwargs):
        """
        Get YambopyBandStructure object with the KS and GW bands
        """
        #create bandstructure objects
        #TODO: should not be kpoints_iku but kpoints_car here
        ks_bandstructure = YambopyBandStructure(self.eigenvalues_dft,self.kpoints_iku,**kwargs)
        qp_bandstructure = YambopyBandStructure(self.eigenvalues_qp, self.kpoints_iku,**kwargs)

        return ks_bandstructure, qp_bandstructure

    def interpolate(self,lattice,path,what='QP+KS',lpratio=5,valence=None,verbose=1,**kwargs):
        """
        Interpolate the QP corrections on a k-point path, requires the lattice structure
        """
        from yambopy.tools.skw import SkwInterpolator

        if verbose:
            print("This interpolation is provided by the SKW interpolator implemented in Abipy")

        cell = (lattice.lat, lattice.red_atomic_positions, lattice.atomic_numbers)
        nelect = 0
        fermie = kwargs.pop('fermie',0)

        #consistency check with lattice from YamboSaveDB
        if len(lattice.kpts_iku)!=len(self.kpoints_iku):
            print(len(lattice.kpts_iku),len(self.kpoints_iku))
            raise ValueError("The QP database is not consistent with the lattice")

        if not np.isclose(lattice.kpts_iku,self.kpoints_iku).all():
            print(lattice.kpts_iku)
            print(self.kpoints_iku)
            raise ValueError("The QP database is not consistent with the lattice")

        #interpolate the dft eigenvalues
        kpoints = lattice.red_kpoints
        sym_rec  = lattice.sym_rec
        symrel = [sym for sym,trev in zip(lattice.sym_rec_red,lattice.time_rev_list) if trev==False ]
        time_rev = True
        
        kpoints_path = path.get_klist()[:,:3]

        #interpolate KS
        ks_ebands, qp_ebands = None, None
        if 'KS' in what:
           if self.spin == True:
              print('Spin-polarized bands DFT')
              eigens_up = self.eigenvalues_dft[np.newaxis,:,:,0]
              eigens_dw = self.eigenvalues_dft[np.newaxis,:,:,1]
              skw_up = SkwInterpolator(lpratio,kpoints,eigens_up,fermie,nelect,cell,symrel,time_rev,verbose=verbose)
              skw_dw = SkwInterpolator(lpratio,kpoints,eigens_dw,fermie,nelect,cell,symrel,time_rev,verbose=verbose)
              dft_eigens_up_kpath = skw_up.interp_kpts(kpoints_path).eigens[0]
              dft_eigens_dw_kpath = skw_dw.interp_kpts(kpoints_path).eigens[0]

              #if valence: kwargs['fermie'] = np.max(dft_eigens_kpath[:,:valence])
              # tricky 
              ks_ebands_up = YambopyBandStructure(dft_eigens_up_kpath,kpoints_path,kpath=path,**kwargs)
              ks_ebands_dw = YambopyBandStructure(dft_eigens_dw_kpath,kpoints_path,kpath=path,**kwargs)

           else:
              print('No spin-polarized bands DFT')
              eigens  = self.eigenvalues_dft[np.newaxis,:]
              skw = SkwInterpolator(lpratio,kpoints,eigens,fermie,nelect,cell,symrel,time_rev,verbose=verbose)
              #kpoints_path = path.get_klist()[:,:3]
              dft_eigens_kpath = skw.interp_kpts(kpoints_path).eigens[0]
              if valence: kwargs['fermie'] = np.max(dft_eigens_kpath[:,:valence])
              ks_ebands = YambopyBandStructure(dft_eigens_kpath,kpoints_path,kpath=path,**kwargs)

        #interpolate QP
        if 'QP' in what:
            if self.spin == True:
               print('Spin-polarized bands QP')
               eigens_up = self.eigenvalues_qp[np.newaxis,:,:,0]
               eigens_dw = self.eigenvalues_qp[np.newaxis,:,:,1]
               print('GW eigenvalues are sorted in ascending energy')
               #sorting
               aux_up, aux_dw = eigens_up, eigens_dw
               for ik in range(self.nkpoints):
                   eigens_up[0,ik,:], eigens_dw[0,ik,:] = sorted(aux_up[0,ik,:]), sorted(aux_dw[0,ik,:])
               #end sorting

               skw_up = SkwInterpolator(lpratio,kpoints,eigens_up,fermie,nelect,cell,symrel,time_rev,verbose=verbose)
               skw_dw = SkwInterpolator(lpratio,kpoints,eigens_dw,fermie,nelect,cell,symrel,time_rev,verbose=verbose)
               #kpoints_path = path.get_klist()[:,:3]
               qp_eigens_up_kpath = skw_up.interp_kpts(kpoints_path).eigens[0]
               qp_eigens_dw_kpath = skw_dw.interp_kpts(kpoints_path).eigens[0]
               #if valence: kwargs['fermie'] = np.max(dft_eigens_kpath[:,:valence])
               # tricky 
               qp_ebands_up = YambopyBandStructure(qp_eigens_up_kpath,kpoints_path,kpath=path,**kwargs)
               qp_ebands_dw = YambopyBandStructure(qp_eigens_dw_kpath,kpoints_path,kpath=path,**kwargs)

            else:
               print('No spin-polarized bands DFT')
               eigens  = self.eigenvalues_qp[np.newaxis,:]
               #sorting
               aux = eigens
               for ik in range(self.nkpoints):
                   eigens[0,ik,:] = sorted(aux[0,ik,:])
               #end sorting
               skw = SkwInterpolator(lpratio,kpoints,eigens,fermie,nelect,cell,symrel,time_rev,verbose=verbose)
               #kpoints_path = path.get_klist()[:,:3]
               qp_eigens_kpath = skw.interp_kpts(kpoints_path).eigens[0]
               if valence: kwargs['fermie'] = np.max(qp_eigens_kpath[:,:valence])

               qp_ebands = YambopyBandStructure(qp_eigens_kpath,kpoints_path,kpath=path,**kwargs)
               #qp_ebands = YambopyBandStructure(qp_eigens_kpath,kpoints_path,kpath=path,weights=qp_z_kpath,size=0.1,**kwargs)

            qp_z_kpath = None
            if 'Z' in what:
                eigens = self.z[np.newaxis,:]
                skw = SkwInterpolator(lpratio,kpoints,eigens,fermie,nelect,cell,symrel,time_rev,verbose=verbose)
                #kpoints_path = path.get_klist()[:,:3]
                qp_z_kpath = skw.interp_kpts(kpoints_path).eigens[0]
                

        if self.spin == True:
           return ks_ebands_up, ks_ebands_dw , qp_ebands_up, qp_ebands_dw
        else: 
           return ks_ebands, qp_ebands

    def expand_eigenvalues(self,lattice):
        """ Expand QP values in full BZ
            - Input: YamboLatticeDB object with expanded kpts
        """
        nkbz   = lattice.nkpoints
        bz2ibz = lattice.kpoints_indexes

        eigenvalues_qp_expanded = np.zeros((nkbz,self.nbands))
        for ikbz in range(nkbz): eigenvalues_qp_expanded[ikbz,:] = self.eigenvalues_qp[bz2ibz[ikbz],:] 
        return eigenvalues_qp_expanded


    @add_fig_kwargs
    def plot_bs(self,**kwargs):
        """
        Get and plot QP bandstructure
        """
        ks_bs,qp_bs = self.get_bs(**kwargs)
        ybs = YambopyBandStructureList([ks_bs,qp_bs])
        return ybs.plot(show=False)

    @property
    def nqps(self):
        return len(self.e)

    @property
    def min_kpoint(self):
        return min(self.kpoint_index)

    @property
    def max_kpoint(self):
        return max(self.kpoint_index)

    @property
    def nbands(self):
        return self.max_band-self.min_band+1

    @property
    def min_band(self):
        return min(self.band_index)

    @property
    def max_band(self):
        return max(self.band_index)

    @property
    def nkpoints(self):
        return len(self.kpoints_iku)
    
    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app("nqps:     %d"%self.nqps)
        app("nkpoints: %d"%self.nkpoints)
        app("nbands:   %d"%self.nbands)
        app("min_band: %d"%self.min_band)
        app("max_band: %d"%self.max_band)
        return "\n".join(lines)

