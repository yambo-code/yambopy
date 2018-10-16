
import numpy as np
from yambopy.plot.plotting import add_fig_kwargs

class YambopyBandStructure():
    """
    Class to plot bandstructures
    """
    _colormap = 'rainbow'

    def __init__(self,bands,kpoints,**kwargs):
        self.bands = bands
        self.kpoints = np.array(kpoints)
        self.kwargs = kwargs
        self._xlim = None
        self._ylim = None

    @property
    def xlim(self):
        if self._xlim is None: return (min(self.distances),max(self.distances))
        return self._xlim

    @property
    def nbands(self):
        nkpoints, nbands = self.bands.shape
        return nbands

    @property
    def nkpoints(self):
        nkpoints, nbands = self.bands.shape
        return nkpoints

    @property
    def ylim(self):
        if self._ylim is None: return (np.min(self.bands),np.max(self.bands))
        return self._ylim

    def set_xlim(self,xlim):
        self._xlim = xlim

    def set_ylim(self,ylim):
        self._ylim = ylim

    @property
    def distances(self):
        if not hasattr(self,"_distances"):
            self._distances = [0]
            distance = 0
            for nk in range(1,len(self.kpoints)):
                distance += np.linalg.norm(self.kpoints[nk]-self.kpoints[nk-1])
                self._distances.append(distance)
        return self._distances

    @add_fig_kwargs    
    def plot(self):
        """return a matplotlib figure with the plot"""
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_ax(ax)
        return fig

    def set_kwargs(self,**kwargs):
        self.kwargs.update(kwargs)

    def get_kwargs(self,**kwargs):
        from copy import deepcopy
        fullkwargs = deepcopy(self.kwargs)
        fullkwargs.update(kwargs)
        return fullkwargs

    def plot_ax(self,ax,xlim=None,ylim=None,legend=False,**kwargs):
        """Receive an intance of matplotlib axes and add the plot"""
        for band in self.bands.T:
            ax.plot(self.distances,band,**self.get_kwargs(**kwargs))
        ax.set_xlim(self.xlim)    
        ax.set_ylim(self.ylim)    
        ax.set_ylabel('Energies (eV)')
        ax.xaxis.set_ticks([])
        if legend: ax.legend()
    
    def __add__(self,y):
        """Add the bands of two systems together"""
        #add some consistency check
        bands = self.bands + y.bands
        return YambopyBandStructure(bands,self.kpoints,**self.kwargs)
 
    def __sub__(self,y):
        """Subtract the bands of two systems together"""
        #add some consistency check
        bands = self.bands - y.bands
        return YambopyBandStructure(bands,self.kpoints,**self.kwargs)

    def __mul__(self,y):
        """Scale the bands of the system"""
        #add some consistency check
        bands = self.bands*y
        return YambopyBandStructure(bands,self.kpoints,**self.kwargs)

    def __truediv__(self,y):
        """Scale the bands of the system"""
        #add some consistency check
        bands = self.bands/y
        return YambopyBandStructure(bands,self.kpoints,**self.kwargs)
    
    def __str__(self):
        lines = []; app = lines.append
        app('nkpoints: %d'%self.nkpoints)
        app('nbands: %d'%self.nbands)
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

    def plot_ax(self,ax,legend=True,**kwargs):
        for bandstructure in self.bandstructures:
            bandstructure.plot_ax(ax)
        if legend and self.has_legend: ax.legend()
 
    @add_fig_kwargs
    def plot(self,**kwargs):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_ax(ax)
        return fig

    def get_colors(self):
        """get a list of colors for each plot"""
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(self._colormap) #get color map
        return [cmap(i) for i in np.linspace(0, 1, len(self.bands))]

    def pickle(self,filename):
        import pickle
        with open(filename,'wb') as f:
            pickle.dump(self,f)

    def __str__(self):
        lines = []; app = lines.append
        app(marqee(self.__cls__.__name__))
        app('nbandstructures: %d'%nbandstructures)
        return "\n".join(lines)
