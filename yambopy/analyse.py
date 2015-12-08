# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
import json
import numpy as np
import re

#we try to use matplotlib, if not present we won't use it
try:
    from matplotlib import pyplot as plt
except ImportError:
    _has_matplotlib = False
else:
    _has_matplotlib = True

class YamboAnalyser():
    """ This class can open multiple yambo files, organize them and plot
    convergence graphs.

    TODO:
    plot convergence of variables according to differences in the input files
    """
    _colormap = 'rainbow'

    def __init__(self,folder='.'):
        self.folder = folder

        files = ["%s/%s"%(folder,filename) for filename in os.listdir(folder)]
        self.filenames = [f for f in files if '.json' in f]

        #read the files
        files = [open(f) for f in self.filenames]
        self.jsonfiles = dict([(filename,json.load(f)) for filename,f in zip(self.filenames,files)])
        for f in files: f.close()

    def get_data_file(self,calculation,tags):
        for filename in self.jsonfiles[calculation]["data"].keys():
            if all(i in filename for i in tags):
                return np.array( self.jsonfiles[calculation]["data"][filename] )

    def get_data(self,tags):
        """ Get a dictionary with all the data from the files under analysis
        """
        data = dict()
        for k in sorted(self.jsonfiles.keys()):
            for filename in self.jsonfiles[k]["data"].keys():
                if all(i in filename for i in tags):
                    data[k] = np.array( self.jsonfiles[k]["data"][filename] )
        return data

    def get_colors(self,tags):
        """ select the colors according to the number of files to plot
        the files to plot are the ones that have all the tags in their name
        """
        nfiles=sum([all(i in filename for i in tags) for k in self.jsonfiles.keys() for filename in self.jsonfiles[k]["data"].keys()])
        cmap = plt.get_cmap(self._colormap) #get color map
        colors = [cmap(i) for i in np.linspace(0, 1, nfiles)]
        return colors

    def plot_gw(self,tags,cols=(lambda x: x[2]+x[3],),rows=None):
        """ Use this function to plot the quasiparticle energies from a GW calculation
            cols: a list of indexes or functions
            the functions that take as input each column of the o.QP file and the output will be plotted
            the index is used to select the columns to plot
            rows: the same as cols but for the electronic bands in the file o.QP

            Example:
                a.plot_gw('qp',cols=(lambda x: x[2]+x[3],),rows=(lambda x: x[1]-x[2],))
                
                Will plot only files with 'qp' in their filename
                Will add the second and third columns (DFT eigenvalues + GW corrections)
                Will subtract the 2nd and 1st bands (usefull to study the convergence of the gap)
        """
        ax = plt.axes([0.1, 0.1, .7, .7])
        plot = False

        colors = self.get_colors(tags)

        n=0
        for k in sorted(self.jsonfiles.keys()):
            for filename in self.jsonfiles[k]["data"].keys():
                if all(i in filename for i in tags):
                    data = np.array( self.jsonfiles[k]["data"][filename] )
                    #the list of quaisparticle energies in yambo is organized as follows:
                    #K-point Band E0 E-E0 Sc(E0)
                    # if we want to plot the bandstructure we have to plot in the same k-point
                    # the values of the required column for the different bands

                    # first we get the number of bands to plot
                    bands = data[:,1]
                    bmin, bmax = int(min(bands)), int(max(bands))
 
                    #select the color to plot with
                    color = colors[n]
                    n+=1
    
                    for col in cols:
                        #get x
                        x = data[data[:,1]==bmin,0]

                        #get the y's
                        #to choose what to plot we can have either a function or a index
                        if hasattr(col, '__call__'):
                            ys = np.array([[ col(c) for c in data[data[:,1]==b,:] ] for b in xrange(bmin,bmax+1)])
                        elif isinstance( col, int ):
                            ys = np.array([ data[data[:,1]==b,col] for b in xrange(bmin,bmax+1) ])
                        else:
                            print "The col datatype: %s is not known"%str(type(col))
                            raise RuntimeError
                       
                        #make row operations if available
                        if rows:
                            ys = [ row(ys) for row in rows ] 

                        #plot 
                        label = filename.split('/')[-1]
                        ax.plot(x,ys[0],'-',label=label,color=color)
                        for y in ys[1:]:
                            ax.plot(x,y,'-',color=color)
                        plot = True
        if plot:
            ax.legend(bbox_to_anchor=(1.05, 1.), loc=2, borderaxespad=0.,prop={'size':8})
            plt.show()

    def plot_bse(self,tags,cols=(2,)):
        """ Use this function to plot the absorption spectrum calculated using the BSE
            cols: a list of indexes to select which columns from the file to plot

            Example:
                a.plot_gw('eps',cols=(2,))
                
                Will plot only files with 'eps' in their filename (absorption spectra)
                Will plot the second column (absorption spectra)
        """
        ax = plt.axes([0.1, 0.1, .7, .7])
        plot = False

        colors = self.get_colors(tags)

        n=0
        for k in sorted(self.jsonfiles.keys()):
            for filename in self.jsonfiles[k]["data"].keys():
                if all(i in filename for i in tags):
                    data = np.array( self.jsonfiles[k]["data"][filename] )

                    #select the color to plot with
                    color = colors[n]
                    n+=1

                    for col in cols:
                        x = data[:,0]
                        y = data[:,col-1]
                        label = filename.split('/')[-1]+"_col=%d"%col
                        ax.plot(x,y,label=label,color=color)
                        plot = True
        if plot:
            ax.legend(bbox_to_anchor=(1.05, 1.), loc=2, borderaxespad=0.,prop={'size':8})
            plt.show()

    def print_timing(self,tags=""):
        for k in self.jsonfiles.keys():
            if all(i in k for i in tags):
                print "\n%s"%k
                for key,val in self.jsonfiles[k]["runtime"].items():
                    print "%40s %10s %10s %10s"%(key,val[0],val[1],val[2])

    def get_inputfiles_tag(self,tags):
        inputfiles = self.get_inputfiles()
        inputfiles_tags = dict()

        for k in inputfiles.keys():
            inputfiles_tags[k] = dict()
            for tag in tags:
                if tag in inputfiles[k].keys():
                    inputfiles_tags[k][tag] = inputfiles[k][tag]
        return inputfiles_tags

    def get_inputfiles(self):
        inputfiles = dict()
       
        for k in self.jsonfiles.keys():
            inputfiles[k] = dict()
            for datatype in ['array','real','string','complex']:
                for key,val in self.jsonfiles[k]["inputfile"][datatype].items():
                    inputfiles[k][key] = val
        return inputfiles

    def print_inputfiles(self):
        for k in self.jsonfiles.keys():
            print "filename:", k
            y = yamboin()
            y.read_variables_dict( self.jsonfiles[k]["inputfile"] )
            print y

    def __str__(self):
        s = ""
        for f in self.filenames:
            s+="%s\n"%f
        return s
