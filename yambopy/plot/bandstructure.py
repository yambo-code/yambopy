
import numpy as np
from yambopy.plot.plotting import add_fig_kwargs

class YamboBandStructure():
    """
    Class to plot bandstructures
    """
    _colormap = 'rainbow'

    def __init__(self,bands=[],distances=[],args=None):
        self.bands = bands
        self.distances = distances

        if args is None:
            self.args = []
        else:
            self.args = args
            

    def add_bands(self,bands,distances=None,**kwargs):
        """
        Add a set of bands to the bandstructure
        
        arguments:
            bands is a one dimensional array with the data to be plotted
        """
        self.bands.append(bands)

        if distances is None: 
            distances = list(range(len(bands)))
        self.distances.append(distances)

        self.args.append(kwargs)

    def get_colors(self):
        """get a list of colors for each plot"""
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(self._colormap) #get color map
        return [cmap(i) for i in np.linspace(0, 1, len(self.bands))]

    @add_fig_kwargs    
    def plot(self):
        """return a matplotlib figure with the plot"""
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_ax(ax)
        return fig        

    def plot_ax(self,ax,xlim=None,ylim=(None,None)):
        """receive an intance of matplotlib axes and add the plot"""
        colors = self.get_colors()
        tmp_xlim = None
        for x,bands,color,args in zip(self.distances,self.bands,colors,self.args):
            for band in bands.T:
                ax.plot(x,band,c=color,**args)
                tmp_xlim = (np.min(x),np.max(x))
                if "label" in args: args.pop("label")
        if xlim is None: xlim = tmp_xlim
        ax.set_xlim(*xlim)    
        if not ylim: ylim = (min(bands),max(bands))
        ax.set_ylim(*ylim)    
        ax.set_ylabel('Energies (eV)')
        ax.xaxis.set_ticks([])
        ax.legend()
    
    def __add__(self,y):
        bands = self.bands + y.bands
        distances = self.distances + y.distances
        args = self.args + y.args
        return YambopyBandStructure(bands=bands,distances=distances,args=args)
    
    def __str__(self):
        s = ""
        s += "number of datasets: %d"%len(self.bands)
        return s
