# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import numpy as np
from yambopy.tools.string import marquee
from yambopy.plot.plotting import add_fig_kwargs
from qepy.lattice import Path

def exagerate_differences(ks_ebandsc,ks_ebandsp,ks_ebandsm,d=0.01,exagerate=5):
    """
    Take three different band-structures with some parameter changing
    and return an exagerated version for plotting.
    This is used for finite differences w.r.t. atomic positions for example
    """
    #calculate central finite difference
    fdc_ks = (ks_ebandsp-ks_ebandsm)/(d*2)
    fdp_ks = (ks_ebandsp-ks_ebandsc)/d
    fdm_ks = (ks_ebandsm-ks_ebandsc)/d

    #exagerate differences
    ks_ebandsp = ks_ebandsc + fdp_ks*d*exagerate
    ks_ebandsm = ks_ebandsc + fdm_ks*d*exagerate

    #set colors again
    ks_ebandsp.set_kwargs(c='brown')
    ks_ebandsm.set_kwargs(c='tomato')

    return ks_ebandsc,ks_ebandsp,ks_ebandsm

def apply_scissor_shift(eigenvalues,scissor,n_val):
    """
    Apply scissor shift to band structure
    
    Input:
      :: eigenvalues -> np array of dimensions (Nk,Nb) or (Ns,Nk,Nb)
      :: scissor -> [shift, stretch_cond, stretch_val]
      :: n_val -> number of valence bands 

      NB: make sure eigenvalues and scissor have same units!        

    Returns shifted eigenvalues
    """
    from copy import deepcopy

    # Dimensionality including spin
    # In this case reshape/concatenate to get (Nkponts,Nbands) array
    if len(eigenvalues.shape)==3:
        Nspin, Nkpoints, Nbands = eigenvalues.shape
        if Nspin>1: raise NotImplementedError("Scissor for spin-polarised bands not yet implemented.")
        else: eigen_to_shift = deepcopy(eigenvalues[0])
     # If original dimensionality is (Nkpoints,Nbands) work with original array
    else: 
        Nkpoints, Nbands = eigenvalues.shape
        eigen_to_shift = deepcopy(eigenvalues)

    # Actual scissor operator code          
    aux = np.zeros((Nkpoints,Nbands))
    top_v, bottom_c = eigen_to_shift[:,n_val-1], eigen_to_shift[:,n_val]
    ind_k_dir_gap = np.argmin(bottom_c-top_v)
    ev_max, ec_min = top_v[ind_k_dir_gap], bottom_c[ind_k_dir_gap]

    for ib in range( Nbands ):
        if ib<n_val: aux[:,ib] = ev_max-(ev_max-eigen_to_shift[:,ib])*scissor[2]
        else:        aux[:,ib] = ec_min+scissor[0]+(eigen_to_shift[:,ib]-ec_min)*scissor[1]

    if len(eigenvalues.shape)==3: final_eigen = np.expand_dims(aux, axis=0)
    else:                         final_eigen = aux

    return final_eigen

class YambopyBandStructure():
    """
    Class to plot bandstructures
    I include spin projection (to be improve and checked)  AMS
    """
    _colormap = 'rainbow'

    def __init__(self,bands,kpoints,kpath=None,fermie=0,weights=None,spin_proj=None,**kwargs):
        self.bands = np.array(bands)
        self.weights = np.array(weights)     if weights   is not None else None
        self.spin_proj = np.array(spin_proj) if spin_proj is not None else None 
        self.kpoints = np.array(kpoints)
        self.kwargs = kwargs
        self.kpath = kpath
        self.fermie = fermie
        self._xlim = None
        self._ylim = None

    @property
    def nbands(self):
        nkpoints, nbands = self.bands.shape
        return nbands

    @property
    def nkpoints(self):
        nkpoints, nbands = self.bands.shape
        return nkpoints

    @property
    def xlim(self):
        if self._xlim is None: return (min(self.distances),max(self.distances))
        return self._xlim

    @property
    def ylim(self):
        if self._ylim is None: return (np.min(self.bands)-self.fermie,np.max(self.bands)-self.fermie)
        return self._ylim

    @classmethod
    def from_dict(cls,d):
        path = Path.from_dict(d['kpath'])
        instance = cls(d['bands'],d['kpoints'],kpath=path,
                       fermie=d['fermie'],weights=d['weights'],**d['kwargs'])
        instance._xlim = d['_xlim']
        instance._ylim = d['_ylim']
        return instance

    @classmethod
    def from_json(cls,filename):
        import json
        with open(filename,'r') as f:
            d = json.load(f)
        return cls.from_dict(d)

    def as_dict(self):
        """ Return the data of this object as a dictionary
        """
        d = { 'bands': self.bands.tolist(),
              'weights': self.weights.tolist() if self.weights is not None else None,
              'kpoints': self.kpoints.tolist(),
              'kwargs': self.kwargs,
              'kpath': self.kpath.as_dict() if self.kpath is not None else None,
              'fermie': self.fermie,
              '_xlim': self._xlim,
              '_ylim': self._ylim }
        return d 

    def write_json(self,filename):
        """serialize this class as a json file"""
        import json
        with open(filename,'w') as f:
            json.dump(self.as_dict(),f)

    def set_fermi(self,valence):
        """simple function to set the fermi energy given the number of valence bands
        """
        self.fermie = np.max(self.bands[:,valence-1])
        self.set_ylim(None)

    def set_energy_offset(self,energy_offset):
        """simple function to rigid-shift the bands
        """
        self.fermie = energy_offset
        self.set_ylim(None)

    def set_xlim(self,xlim):
        self._xlim = xlim

    def set_ylim(self,ylim):
        self._ylim = ylim

    def set_ax_lim(self,ax,fermie=0,ylim=None,xlim=None):
        if xlim is None: xlim = self.xlim
        if ylim is None: ylim = self.ylim
        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0]-self.fermie,ylim[1]-self.fermie)

    @property
    def distances(self):
        if not hasattr(self,"_distances"):
            self._distances = [0]
            distance = 0
            for nk in range(1,len(self.kpoints)):
                distance += np.linalg.norm(self.kpoints[nk]-self.kpoints[nk-1])
                self._distances.append(distance)
        return self._distances

    def as_list(self,bands=None):
        yl = YambopyBandStructureList([self])
        if bands: yl.append(bands)
        return yl

    @add_fig_kwargs    
    def plot(self,title=None):
        """return a matplotlib figure with the plot"""
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_ax(ax)
        if title: ax.title(title)
        return fig

    def set_kwargs(self,**kwargs):
        self.kwargs.update(kwargs)

    def get_kwargs(self,**kwargs):
        from copy import deepcopy
        fullkwargs = deepcopy(self.kwargs)
        fullkwargs.update(kwargs)
        return fullkwargs

    def add_kpath_labels(self,ax):
        """
        Add vertical lines at the positions of the high-symmetry k-points
        I did a modification that I don't like much. To be corrected
        """
        if self.kpath is None:
            ax.xaxis.set_ticks([])
            return 
        for kpoint, klabel, distance in self.kpath:
            ax.axvline(distance,c='k',ls='--',lw=0.5)
        ax.axvline(0.0,c='k',ls='-',lw=0.0)
        ax.axvline(distance,c='k',ls='-',lw=0.0)
        self.kpath.set_xticks(ax)

    def plot_ax(self,ax,xlim=None,ylim=None,size=1.,ylabel='$\epsilon_{n\mathbf{k}}$ [eV]', alpha_weights=0.5,legend=False,**kwargs):
        """Receive an intance of matplotlib axes and add the plot"""
        import matplotlib.pyplot as plt
        kwargs = self.get_kwargs(**kwargs)
        fermie = kwargs.pop('fermie',self.fermie)

        # Set kwargs
        c_bands   = kwargs.pop('c_bands',None)
        c_weights = kwargs.pop('c_weights',None)
        label   = kwargs.pop('label',None)
        lw_label  = kwargs.pop('lw_label',None)
        marker = kwargs.pop('marker',None)
        linestyle = kwargs.pop('linestyle',None)

        # I choose a colormap for spin
        color_map  = plt.get_cmap('seismic')

        x = self.distances
        y = self.bands.T-fermie
        for ib,band in enumerate(self.bands.T):
            x = self.distances
            y = band-fermie
            ax.plot(x,y,c=c_bands,lw=lw_label,marker=marker,linestyle=linestyle,label=label if ib == 0 else "_nolegend_")
            # fill between 
            if self.weights is not None:

                dy = self.weights[:,ib]*size
                ax.fill_between(x,y+dy,y-dy,alpha=alpha_weights,color=c_weights,linewidth=0,label=label)

        self.set_ax_lim(ax,fermie=fermie,xlim=xlim,ylim=ylim)
        ax.set_ylabel(ylabel)
        self.add_kpath_labels(ax)
        if legend: ax.legend()

    def plot_spin_ax(self,ax,xlim=None,ylim=None,ylabel='$\epsilon_{n\mathbf{k}}$[eV]',alpha_weights=0.5,spin_proj_bands=None,legend=False,**kwargs):
        """Receive an intance of matplotlib axes and add the plot"""
        #
        # There is a problem with the number of points in the k-poitns path
        import matplotlib.pyplot as plt
        kwargs = self.get_kwargs(**kwargs)
        fermie = kwargs.pop('fermie',self.fermie)

        # Set kwargs
        c_bands   = kwargs.pop('c_bands',None)
        c_weights = kwargs.pop('c_weights',None)
        label   = kwargs.pop('label',None)
        lw_label  = kwargs.pop('lw_label',None)
        marker = kwargs.pop('marker',None)
        linestyle = kwargs.pop('linestyle',None)
        n_valence = kwargs.pop('n_valence',None)
        weight_option = kwargs.pop('weight_option',None)

        # I choose a colormap for spin
        color_map  = plt.get_cmap('PiYG')

        x = self.distances
        y = self.bands.T-fermie
        nk_distance = len(self.distances)
        for ib,band in enumerate(self.bands.T):
            x = self.distances
            y = band-fermie
            ax.plot(x,y,c=c_bands,lw=lw_label,marker=marker,linestyle=linestyle,label=label if ib == 0 else "_nolegend_")

            # dot
            if self.weights is not None:
               dy = self.weights[:,ib]*size*1000
               color_spin = spin_proj_bands[:,n_valence + ib] + 0.5 # I renormalize 0 => down; 1 => up
               ax.scatter(x,y,s=abs(dy),c=color_spin,cmap=color_map,edgecolors='none',zorder=2,rasterized=True) 

        self.set_ax_lim(ax,fermie=fermie,xlim=xlim,ylim=ylim)
        ax.set_ylabel(ylabel)
        self.add_kpath_labels(ax)
        if legend: ax.legend()

    def __add__(self,y):
        """Add the bands of two systems together"""
        #add some consistency check
        bands = self.bands + y.bands
        return YambopyBandStructure(bands,self.kpoints,kpath=self.kpath,fermie=self.fermie+y.fermie,**self.kwargs)
 
    def __sub__(self,y):
        """Subtract the bands of two systems together"""
        #add some consistency check
        bands = self.bands - y.bands
        return YambopyBandStructure(bands,self.kpoints,kpath=self.kpath,fermie=self.fermie-y.fermie,**self.kwargs)

    def __mul__(self,y):
        """Scale the bands of the system"""
        #add some consistency check
        bands = self.bands*y
        return YambopyBandStructure(bands,self.kpoints,kpath=self.kpath,fermie=self.fermie*y,**self.kwargs)

    def __truediv__(self,y):
        """Scale the bands of the system"""
        #add some consistency check
        bands = self.bands/y
        return YambopyBandStructure(bands,self.kpoints,kpath=self.kpath,fermie=self.fermie/y,**self.kwargs)
    
    def __str__(self):
        lines = []; app = lines.append
        app('nkpoints: %d'%self.nkpoints)
        app('nbands: %d'%self.nbands)
        app('has kpath: %s'%hasattr(self,'kpath'))
        return "\n".join(lines)

class YambopyBandStructureList():
    """This class contains a list of band-structure classes and is responsible for plotting them together"""
    def __init__(self,bandstructures):
        self.bandstructures = bandstructures

    @classmethod
    def from_pickle(cls,filename):
        import pickle
        with open(filename,'rb') as f:
            ybs = pickle.load(f)
        return ybs

    @property
    def has_legend(self):
        return any(['label' in bandstructure.kwargs for bandstructure in self.bandstructures])

    @property
    def nbandstructures(self):
        return len(self.bandstructures)

    @property
    def xlim(self):
        low_xlim = [bandstructure.xlim[0] for bandstructure in self.bandstructures]
        top_xlim = [bandstructure.xlim[1] for bandstructure in self.bandstructures]
        return (np.min(low_xlim),np.max(top_xlim))

    @property
    def ylim(self):
        low_ylim = [bandstructure.ylim[0] for bandstructure in self.bandstructures]
        top_ylim = [bandstructure.ylim[1] for bandstructure in self.bandstructures]
        return (np.min(low_ylim),np.max(top_ylim))

    def as_dict(self):
        bandstructures_list=[]
        for bandstructure in self.bandstructures:
            bandstructures_list.append(bandstructure.as_dict())
        return bandstructures_list

    @classmethod
    def from_json(cls,filename):
        import json
        with open(filename,'r') as f:
            d = json.load(f)
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls,d):
        bandstructures = []
        for bandstructure_dict in d:
            bandstructure = YambopyBandStructure.from_dict(bandstructure_dict)
            bandstructures.append(bandstructure)
        return cls(bandstructures)

    def __getitem__(self,idx):
        return self.bandstructures[idx]
    
    def append(self,bands):
        if isinstance(bands,list): self.bandstructure.extend(bands)
        self.bandstructures.append(bands)

    def add_bandstructure(self,bandstructures,**kwargs):
        """
        Add a bandstructure to bandstructure set
        """
        if not isinstance(bandstructures,list): bandstructures = [bandstructures]

        #add arguments in all the structures
        for bandstructure in bandstructures:
            bandstructure.kwargs.update(kwargs)

        #extend current bandstructure set
        self.bandstructures.extend(bandstructures)

    def plot_ax(self,ax,legend=True,xlim=None,ylim=None,**kwargs):
        title = kwargs.pop('title',None)
        for i,bandstructure in enumerate(self.bandstructures):
            bandstructure.plot_ax(ax,**kwargs)
        if xlim is None: ax.set_xlim(self.xlim)
        else:            ax.set_xlim(xlim)
        if ylim is None: ax.set_ylim(self.ylim)
        else:            ax.set_ylim(ylim)
        if title: ax.set_title(title)
        if legend and self.has_legend: ax.legend()

    def set_fermi(self,valence):
        """Find the Fermi energy of all the bandstructures by shifting the Fermi energy"""
        for bandstructure in self.bandstructures:
            bandstructure.set_fermi(valence)

    @add_fig_kwargs
    def plot(self,figsize=None,**kwargs):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1,1,1)
        self.plot_ax(ax,**kwargs)
        return fig
    
    def get_color(self,i,colormap='gist_rainbow'):
        colormap = self.get_colormap(colormap=colormap) 
        return colormap[i]

    def get_colormap(self,colormap='gist_rainbow'):
        """get a list of colors for each plot"""
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(colormap) #get color map
        return [cmap(i) for i in np.linspace(0, 1, self.nbandstructures)]

    def write_json(self,filename):
        """serialize this class as a json file"""
        import json
        with open(filename,'w') as f:
            json.dump(self.as_dict(),f)

    def pickle(self,filename):
        import pickle
        with open(filename,'wb') as f:
            pickle.dump(self,f)

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app('nbandstructures: %d'%self.nbandstructures)
        return "\n".join(lines)
