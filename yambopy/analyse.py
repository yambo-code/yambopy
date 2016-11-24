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
from itertools import product
from yambopy import *

#we try to use matplotlib, if not present we won't use it
try:
    from matplotlib import pyplot as plt
except ImportError:
    _has_matplotlib = False
else:
    _has_matplotlib = True

def red_car(red,lat):
    lat = np.array(lat)
    return np.array(map( lambda coord: coord[0]*lat[0]+coord[1]*lat[1]+coord[2]*lat[2], red))

def car_red(car,lat): return np.array(map( lambda coord: np.linalg.solve(np.array(lat).T,coord), car))
def rec_lat(lat):
    """
    Calculate the reciprocal lattice vectors
    """
    a1,a2,a3 = np.array(lat)
    v = np.dot(a1,np.cross(a2,a3))
    b1 = np.cross(a2,a3)/v
    b2 = np.cross(a3,a1)/v
    b3 = np.cross(a1,a2)/v
    return np.array([b1,b2,b3])

def expand_kpts(kpts,syms):
    """ Take a list of qpoints and symmetry operations and return the full brillouin zone
    with the corresponding index in the irreducible brillouin zone
    """
    full_kpts = []
    print "nkpoints:", len(kpts)
    for nk,k in enumerate(kpts):
        for sym in syms:
            full_kpts.append((nk,np.dot(sym,k)))

    return full_kpts

def isbetween(a,b,c):
    """ Check if c is between a and b
    """
    return np.isclose(np.linalg.norm(a-c)+np.linalg.norm(b-c)-np.linalg.norm(a-b),0)

class YamboAnalyser():
    """ This class can open multiple yambo files, organize them and plot
    convergence graphs.

    TODO:
    plot convergence of variables according to differences in the input files
    """
    _colormap = 'rainbow'

    def __init__(self,folder='.'):
        self.folder = folder

        files = ["%s"%filename for filename in os.listdir(folder)]
        self.filenames = [f for f in files if '.json' in f]

        #read the files
        files = [open("%s/%s"%(self.folder,f)) for f in self.filenames]
        self.jsonfiles = dict([(filename,json.load(f)) for filename,f in zip(self.filenames,files)])
        for f in files: f.close()

    def get_data_file(self,calculation,tags):
        for filename in self.jsonfiles[calculation]["data"].keys():
            if all(i in filename for i in tags):
                tags = self.jsonfiles[calculation]["tags"][filename]
                data = np.array( self.jsonfiles[calculation]["data"][filename] )
                return dict( zip(tags, data.T) )

    def get_data(self,tags):
        """ Get a dictionary with all the data from the files under analysis
        """
        data = dict()
        for k in sorted(self.jsonfiles.keys()):
            for filename in self.jsonfiles[k]["data"].keys():
                if all(i in filename for i in tags):
                    data[k] = np.array( self.jsonfiles[k]["data"][filename] )
        return data

    def get_tags(self,tags):
        """ Get a dictionary with the tags of the output file colomns
        """
        tagslist = dict()
        for k in sorted(self.jsonfiles.keys()):
            for filename in self.jsonfiles[k]["tags"].keys():
                if all(i in filename for i in tags):
                    tagslist[k] = np.array( self.jsonfiles[k]["tags"][filename] )
        return tagslist

    def get_colors(self,tags):
        """ select the colors according to the number of files to plot
        the files to plot are the ones that have all the tags in their name
        """
        nfiles=sum([all(i in filename for i in tags) for k in self.jsonfiles.keys() for filename in self.jsonfiles[k]["data"].keys()])
        cmap = plt.get_cmap(self._colormap) #get color map
        colors = [cmap(i) for i in np.linspace(0, 1, nfiles)]
        return colors

    def get_path(self,path,json_filename):
        """ Obtain a list of indexes and kpoints that belong to the regular mesh
        """
        jsonfile = self.jsonfiles[json_filename]

        if 'kpts_iku' not in jsonfile or 'sym_car' not in jsonfile:
            raise ValueError( "Could not find information about the k points in the json file %s."%json_filename )

        #get data from json file
        kpts_iku = np.array(jsonfile['kpts_iku'])
        sym_car  = np.array(jsonfile['sym_car'])
        alat     = np.array(jsonfile['alat'])
        lattice  = np.array(jsonfile['lattice'])

        #check if the lattice data is present
        if not lattice.any():
            raise ValueError('Information about the lattice is not present, cannot determine the path')

        #convert to cartesian coordinates
        kpts_car = np.array([ k/alat for k in kpts_iku ])

        #get the full list of kpoints
        full_kpts = expand_kpts(kpts_car,sym_car)

        #points in cartesian coordinates
        reciprocal_lattice = rec_lat(lattice)
        path_car = red_car(path, reciprocal_lattice)

        #find the points along the high symmetry lines
        distance = 0
        bands_kpoints = []
        bands_indexes = []
        bands_highsym_qpts = []
        old_kpt = np.array([0,0,0])
        for k in range(len(path)-1):

            kpoints_in_path = {} #store here all the points in the path
            start_kpt = path_car[k]
            end_kpt = path_car[k+1]
            #find the collinear points
            for x,y,z in product(range(-1,2),repeat=3):
                shift = red_car(np.array([[x,y,z]]),reciprocal_lattice)[0]
                for index, kpt in full_kpts:
                    kpt_shift = kpt+shift
                    #if the point is collinear and not in the list we add it
                    if isbetween(start_kpt,end_kpt,kpt_shift):
                        key = tuple([round(kpt,4) for kpt in kpt_shift])
                        value = [ index, np.linalg.norm(start_kpt-kpt_shift), kpt_shift ]
                        kpoints_in_path[key] = value

            #sort the points acoording to distance to the start of the path
            kpoints_in_path = sorted(kpoints_in_path.values(),key=lambda i: i[1])

            #get kpoints_in_pathpoints
            if k==0: bands_highsym_qpts.append(kpoints_in_path[0][2])
            for index, disp, kpt in kpoints_in_path:
                bands_kpoints.append( kpt )
                bands_indexes.append( index )
                print ("%12.8lf "*3)%tuple(kpt), index
            bands_highsym_qpts.append(kpt)
        return bands_kpoints, bands_indexes, bands_highsym_qpts

    def plot_gw_path(self,tags,path_label,cols=(lambda x: x[2]+x[3],),rows=None):
        """ Create a path of k-points and find the points in the regular mesh that correspond to points in the path
            Use these points to plot the GW band structure.
        """
        path = np.array([p[0] for p in path_label])
        labels = [p[1] for p in path_label]
        plot = False
        colors = self.get_colors(tags)
        lstyles = ['-', '--', '_', ':']
        fig = plt.figure()
        ax = plt.subplot(111)
        n=0

        #select one of the files to obtain the points in the path
        json_filename = self.jsonfiles.keys()[0]

        #find the points along the high symmetry lines
        json_filename = self.jsonfiles.keys()[0]
        bands_kpoints, bands_indexes, bands_highsym_qpts = self.get_path(path,json_filename)

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk-1]-bands_kpoints[nk])
            bands_distances.append(distance)
        

        #obtain the bands for the output files and plot
        for json_filename in self.jsonfiles.keys():
            for output_filename in self.jsonfiles[json_filename]['data']:
                kpoint_index, bands_cols = self.get_gw_bands(json_filename,output_filename,cols=cols,rows=rows)

                #plot
                for ib,bands in enumerate(bands_cols):
                    label = output_filename
                    for band in bands:
                        plt.plot(bands_distances,[band[k] for k in bands_indexes],linestyle=lstyles[ib%len(lstyles)],label=label,color=colors[n])
                        label=None
                plot = True
                n+=1

        if plot:
            #plot high-symmetry q-points
            distance = 0
            bands_highsym_qpts_distances = [0]
            for nk in range(1,len(bands_highsym_qpts)):
                plt.axvline(distance,color='k')
                distance+=np.linalg.norm(bands_highsym_qpts[nk]-bands_highsym_qpts[nk-1])
                bands_highsym_qpts_distances.append(distance)

            #plot labels
            plt.xticks(bands_highsym_qpts_distances, labels)

            box = ax.get_position()
            plt.title('GW quasiparticles on a path')
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            plt.xlim(0,max(bands_distances))
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
            plt.show()

    def get_gw_bands(self,json_filename,output_filename,cols=(lambda x: x[2]+x[3],),rows=None):
        """ Get the gw bands from a gw calculation from a filename
            json_filename the name of the json file
            output_filename the name of the output filename that is in the json file

            The list of quaisparticle energies in yambo is organized as follows:
            K-point Band E0 E-E0 Sc(E0)
            if we want to plot the bandstructure we have to plot in the same k-point
            the values of the required column for the different bands
        """
        data = np.array( self.jsonfiles[json_filename]["data"][output_filename] )
        # first we get the number of bands to plot
        bands = data[:,1]
        bmin, bmax = int(min(bands)), int(max(bands))

        bands_cols = []
        for col in cols:
            #get x
            kpoint_index = data[data[:,1]==bmin,0]

            #get the y's
            #to choose what to plot we can have either a function or an index
            if hasattr(col, '__call__'):
                bands = np.array([[ col(c) for c in data[data[:,1]==b,:] ] for b in xrange(bmin,bmax+1)])
            elif isinstance( col, int ):
                bands = np.array([ data[data[:,1]==b,col] for b in xrange(bmin,bmax+1) ])
            else:
                raise ValueError( "The col datatype: %s is not known"%str(type(col)) )

            #make row operations if available
            if rows:
                bands = [ row(bands) for row in rows ]
            bands_cols.append(bands)
        return kpoint_index, bands_cols

    def plot_qp_correction(self,tags=('qp',),lda=2,qp=3):
       ax = plt.axes([0.1, 0.1, .7, .7])
       for json_filename in sorted(self.jsonfiles.keys()):
            for output_filename in self.jsonfiles[json_filename]["data"]:
                if all(i in output_filename for i in tags):
                    data = np.array( self.jsonfiles[json_filename]["data"][output_filename] )

                    plt.plot(data[:,lda],data[:,qp],'o',label=output_filename)
                    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
       xmin, xmax = ax.get_xlim()
       ymin, ymax = ax.get_ylim()
       plt.plot([xmin,xmax],[ymin,ymax],'k--',lw=2)

       box = ax.get_position()
       ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
       ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
       plt.show()
       plt.plot()
       plt.show()

    def plot_gw(self,tags=('qp',),cols=(lambda x: x[2]+x[3],),rows=None):
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
        plot = False
        fig = plt.figure()
        ax = plt.subplot(111)
        colors = self.get_colors(tags)

        n=0
        for json_filename in sorted(self.jsonfiles.keys()):
            for output_filename in self.jsonfiles[json_filename]["data"]:
                if all(i in output_filename for i in tags):
                    data = np.array( self.jsonfiles[json_filename]["data"][output_filename] )

                    kpoint_index, bands_cols = self.get_gw_bands(json_filename,output_filename,cols=cols,rows=rows)

                    #plot
                    for bands in bands_cols:
                        label = output_filename
                        for band in bands:
                            ax.plot(kpoint_index,band,'-',label=label,color=colors[n])
                            label = None
                    plot = True
                    n+=1
        if plot:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

            plt.title('GW quasiparticles on a mesh')
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
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
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
            plt.show()

    def plot_spectral_function(self,tags):
        if type(tags) == str:
            tags = (tags,)
        ax = plt.axes([0.1, 0.1, .7, .7])
        for json_filename in sorted(self.jsonfiles.keys()):
          for output_filename in self.jsonfiles[json_filename]["data"]:
            if all(i in output_filename for i in tags):
               data = np.array( self.jsonfiles[json_filename]["data"][output_filename] )
               plt.title('Spectral function as a function of Temperature')
               plt.plot(data[:,0],data[:,2],'-',label=output_filename)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
        plt.show()

    def print_timing(self,tags=""):
        for k in self.jsonfiles.keys():
            if all(i in k for i in tags):
                print "\n%s"%k
                for key,val in self.jsonfiles[k]["runtime"].items():
                    print "%40s %10s %10s %10s"%(key,val[0],val[1],val[2])

    def get_inputfiles_tag(self,tags):
        """
        Get a specific tag from all the .json files from the folders
        You need to write down all the tags that you want to find
        The tags are both for variables in the input file and arguments (meaning runlevels)
        """
        #check if a string was passed and in that case we make it a tuple
        if type(tags) == str:
            tags = (tags,)

        inputfiles = self.get_inputfiles()
        inputfiles_tags = dict()

        for k in inputfiles.keys():
            inputfiles_tags[k] = dict()
            
            # get the current inputfile
            this_inputfile = inputfiles[k]
            
            #initialize the dictionary
            inputfiles_tags[k] = {'variables':{},'arguments':[]}

            for tag in tags:
                for filename in this_inputfile:

                    # We look for the tag both in the variable and in the arguments
                    # look in variables
                    if tag in this_inputfile[filename]['variables'].keys():
                        inputfiles_tags[k]['variables'][tag] = this_inputfile[filename]['variables'][tag]
                    #look in arguments
                    if tag in this_inputfile[filename]['arguments']:
                        inputfiles_tags[k]['arguments'].append(tag)
        return inputfiles_tags

    def get_inputfiles(self):
        """
        Get all the inputfiles from the different .json files
        Each .json file contains all the output files in a folder
        """

        inputfiles = dict()
        for k in self.jsonfiles:
            inputfiles[k] = dict()
            for key,val in self.jsonfiles[k]["inputfile"].items():
                inputfiles[k][key] = val
        return inputfiles

    def print_inputfiles(self):
        """
        Print all the inputfiles from all the json files
        """
        #iterate over the json files
        for k in self.jsonfiles.keys():
            print "jsonfile: ", k
            
            #get the jsonfile
            jsonfile = self.jsonfiles[k]

            for inputfile,content in jsonfile['inputfile'].items():
                print "inputfile:", inputfile
                y = YamboIn(filename=None)
                y.arguments = content["arguments"]
                y.variables = content["variables"]
                print y
                print

    def __str__(self):
        s = ""
        for json_file in self.jsonfiles:
            s+="%s\n"%json_file
            for f in self.jsonfiles[json_file]['data']:
                s+="\t%s\n"%f
        return s
