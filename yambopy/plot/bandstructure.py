
import numpy as np
from yambopy.plot.plotting import add_fig_kwargs

class YamboBandStructure():
    """
    Class to plot bandstructures
    """
    _colormap = 'rainbow'

    def __init__(self,bands=None,distances=None,args=None):
        if bands is None: 
            self.bands = []
        else:
            self.bands = bands

        if distances is None:
            self.distances = []            
        else:
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
            distances = range(len(bands))
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
        for x,bands,color,args in zip(self.distances,self.bands,colors,self.args):
            for band in bands.T:
                ax.plot(x,band,c=color,**args)
                if "label" in args: args.pop("label")
        if not xlim: xlim = (min(x),max(x))
        ax.set_xlim(*xlim)    
        if not ylim: ylim = (min(bands),max(bands))
        ax.set_ylim(*ylim)    
        ax.legend()
    
    def save_pdf(self,filename):
        """save the plot as a pdf file"""
        fig = self.plot()
        plt.savefig("%s.pdf",filename)

    def __add__(self,y):
        bands = self.bands + y.bands
        distances = self.distances + y.distances
        args = self.args + y.args
        return YambopyBandStructure(bands=bands,distances=distances,args=args)
    
    def __str__(self):
        s = ""
        s += "number of datasets: %d"%len(self.bands)
        return s
