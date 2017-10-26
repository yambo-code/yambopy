from __future__ import print_function, division
#
# Copyright (C) 2017 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
import json
import re
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
from yambopy.lattice import red_car, rec_lat, expand_kpts, isbetween
from yambopy.io.inputfile import YamboIn
from yambopy.duck import isstring

class YamboAnalyser():
    """
    Class to open multiple ``.json`` files, organize them and plot the data together.
    Used to perform convergence tests
    """
    _colormap = 'rainbow'

    def __init__(self, folder='.'):
        self.folder = folder

        #get all the json files in the folder
        all_filenames  = os.listdir(folder)

        #filter for json files
        json_filenames = [f for f in all_filenames if '.json' in f]

        #read all the json files
        self.jsonfiles = dict()
        for json_filename in json_filenames:
            path = "%s/%s"%(folder, json_filename)
            
            with  open(path,'r') as f:
                self.jsonfiles[json_filename] = json.load(f)

    def _getter(self,tags,what):
        """
        get data from the json files 
        """
        if isstring(tags):
            tags=(tags,)

        data = dict()
        for k in sorted(self.jsonfiles.keys()):
            for filename in list(self.jsonfiles[k][what].keys()):
                if any(i in filename for i in tags):
                    data[k] = np.array( self.jsonfiles[k][what][filename] )
        return data

    def get_data(self,tags):
        """ 
        Get a dictionary with all the data from the files under analysis
        """
        return self._getter(tags,"data")

    def get_tags(self,tags):
        """ 
        Get a dictionary with the tags of the output file colomns
        """
        return self._getter(tags,"tags")

    def get_colors(self,tags):
        """ 
        Select the colors according to the number of files to plot
        the files to plot are the ones that have all the tags in their name
        """
        #count the number of files
        nfiles = 0
        for k in self.jsonfiles.keys():
            for filename in list(self.jsonfiles[k]["data"].keys()):
                nfiles+=all(i in filename for i in tags)

        cmap = plt.get_cmap(self._colormap) #get color map
        colors = [cmap(i) for i in np.linspace(0, 1, nfiles)]
        return colors

    def get_inputfiles_tag(self,tags):
        """
        Get a specific tag from all the .json files from the folders
        You need to write down all the tags that you want to find
        The tags are both for variables in the input file and arguments (runlevels)
        """
        #check if a string was passed and in that case we make it a tuple
        if isinstance(tags,str):
            tags = (tags,)

        inputfiles = self.get_inputfiles()
        inputfiles_tags = dict()

        for k in list(inputfiles.keys()):
            inputfiles_tags[k] = dict()

            # get the current inputfile
            this_inputfile = inputfiles[k]

            #initialize the dictionary
            inputfiles_tags[k] = {'variables':{},'arguments':[]}

            for tag in tags:
                for filename in this_inputfile:
                    
                    # look in variables
                    if tag in list(this_inputfile[filename]['variables'].keys()):
                        variables = this_inputfile[filename]['variables'][tag]
                        inputfiles_tags[k]['variables'][tag] = variables

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
            for key,val in list(self.jsonfiles[k]["inputfile"].items()):
                inputfiles[k][key] = val
        return inputfiles

    def get_path(self,path,json_filename):
        """
        Obtain a list of indexes and kpoints that belong to the regular mesh
        """
        jsonfile = self.jsonfiles[json_filename]

        if 'kpts_iku' not in jsonfile or 'sym_car' not in jsonfile:
            raise ValueError( "Could not find information about the "
                              "k points in the json file %s."%json_filename )

        #get data from json file
        kpts_iku = np.array(jsonfile['kpts_iku'])
        sym_car  = np.array(jsonfile['sym_car'])
        alat     = np.array(jsonfile['alat'])
        lattice  = np.array(jsonfile['lattice'])

        #check if the lattice data is present
        if not lattice.any():
            raise ValueError('Information about the lattice is not present, '
                             'cannot determine the path')

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
            for x,y,z in product(list(range(-1,2)),repeat=3):
                shift = red_car(np.array([[x,y,z]]),reciprocal_lattice)[0]
                for index, kpt in full_kpts:
                    kpt_shift = kpt+shift
                    #if the point is collinear and not in the list we add it
                    if isbetween(start_kpt,end_kpt,kpt_shift):
                        key = tuple([round(kpt,4) for kpt in kpt_shift])
                        value = [ index, np.linalg.norm(start_kpt-kpt_shift), kpt_shift ]
                        kpoints_in_path[key] = value

            #sort the points acoording to distance to the start of the path
            kpoints_in_path = sorted(list(kpoints_in_path.values()),key=lambda i: i[1])

            #get kpoints_in_pathpoints
            if k==0: bands_highsym_qpts.append(kpoints_in_path[0][2])
            for index, disp, kpt in kpoints_in_path:
                bands_kpoints.append( kpt )
                bands_indexes.append( index )
                print(("%12.8lf "*3)%tuple(kpt), index)
            bands_highsym_qpts.append(kpt)
        return bands_kpoints, bands_indexes, bands_highsym_qpts

    def plot_gw_path(self,path_label,tags=('qp',),type_calc=('lda',),set_calc=(),rows=None):
        """
        Create a path of k-points and find the points in the 
        regular mesh that correspond to points in the path
        Use these points to plot the GW band structure.

        """
        '''
        Assign type of line
        '''
        line_dict  = {'lda':'-','gw':'--','corr':'-.'}
        '''
        Assign color with the number of calculations
        '''
        n_calculations = len(self.jsonfiles.keys())
        cmap = plt.get_cmap(self._colormap) #get color map
        colors = [cmap(i) for i in np.linspace(0, 1, n_calculations)]

        path = np.array([p[0] for p in path_label])
        labels = [p[1] for p in path_label]
        plot = False

        fig = plt.figure()
        ax = plt.subplot(111)
        n=0

        #find the points along the high symmetry lines
        json_filename = list(self.jsonfiles.keys())[0]
        bands_kpoints, bands_indexes, bands_highsym_qpts = self.get_path(path,json_filename)

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk-1]-bands_kpoints[nk])
            bands_distances.append(distance)

        # If we select a set of data
        if set_calc:
          for set_plot in set_calc:
            for i_json,json_filename in enumerate(self.jsonfiles.keys()):
              if set_plot in json_filename:
                for i_type,e_data in enumerate(type_calc):
                  kpoint_index, bands_cols = self.get_gw_bands(json_filename,e_data,json_filename[:-5])
                  label_data = '%s-%s' % (json_filename[:-5],e_data)
                  for ib,bands in enumerate(bands_cols):
                    for i_band,band in enumerate(bands):
                      if i_band == 0:
                        plt.plot(bands_distances,[band[k] for k in bands_indexes],linestyle=line_dict[e_data],label=label_data,color=colors[i_json])
                      else:
                        plt.plot(bands_distances,[band[k] for k in bands_indexes],linestyle=line_dict[e_data],color=colors[i_json])
                      plot = True
                      n+=1

        # If we plot all data
        if not set_calc:
          for i_json,json_filename in enumerate(self.jsonfiles.keys()):
            for i_type,e_data in enumerate(type_calc):
                kpoint_index, bands_cols = self.get_gw_bands(json_filename,e_data,json_filename[:-5])
                label_data = '%s-%s' % (json_filename[:-5],e_data)
                for ib,bands in enumerate(bands_cols):
                  for i_band,band in enumerate(bands):
                    if i_band == 0:
                      plt.plot(bands_distances,[band[k] for k in bands_indexes],linestyle=line_dict[e_data],label=label_data,color=colors[i_json])
                  else:
                      plt.plot(bands_distances,[band[k] for k in bands_indexes],linestyle=line_dict[e_data],color=colors[i_json])
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
            plt.ylabel('E (eV)')
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            plt.xlim(0,max(bands_distances))
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
            if 'DISPLAY' in os.environ:
                plt.show()
            else:
                plt.savefig('gw.png')

    def get_gw_path_bands(self,path_label,tags=('qp',),type_calc=('lda',),rows=None):

        """ Get the bands a path of k-points and find the points in the regular mesh that correspond to points in the path
        """
        path = np.array([p[0] for p in path_label])

        #find the points along the high symmetry lines
        json_filename = list(self.jsonfiles.keys())[0]
        bands_kpoints, bands_indexes, bands_highsym_qpts = self.get_path(path,json_filename)

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk-1]-bands_kpoints[nk])
            bands_distances.append(distance)

        #obtain the bands for the output files and plot
        for json_filename in self.jsonfiles.keys():
            #for output_filename in self.jsonfiles[json_filename]['data']:
            for e_data in type_calc:
                kpoint_index, bands_cols = self.get_gw_bands(json_filename,e_data,json_filename[:-5])
                # Pass path and bands in arrays
                band_in_path = []
                for ib,bands in enumerate(bands_cols):
                    band_aux = []
                    for band in bands:
                        band_aux.append([band[k] for k in bands_indexes])
                    band_in_path.append(band_aux)

        return bands_distances, band_in_path

    def get_gw_bands(self,json_filename,type_calc,output_filename):
        """
        Get the gw bands from a gw calculation from a filename

        Arguments:
            json_filename: the name of the json file
            output_filename: the name of the output filename that is in the json file
            type_calc:
               We read from the netcdf file:
               Eo : LDA
               E  : GW
               E-Eo : GW corrections
        """
        data = np.array( self.jsonfiles[json_filename]["data"][output_filename] )
        b_index = self.jsonfiles[json_filename]["tags"][output_filename].index('Band')
        k_index = self.jsonfiles[json_filename]["tags"][output_filename].index('Kpoint_index')
        Eo_index = self.jsonfiles[json_filename]["tags"][output_filename].index('Eo')
        E_index  = self.jsonfiles[json_filename]["tags"][output_filename].index('E')
        DE_index = self.jsonfiles[json_filename]["tags"][output_filename].index('E-Eo')
        # first we get the number of bands to plot
        bands = data[:,b_index] # new format
        bmin, bmax = int(min(bands)), int(max(bands))

        bands_cols = []

        kpoint_index = data[data[:,b_index]==bmin,k_index]

        #if 'lda' in evalues:
        if type_calc == 'lda':
          bands = [data[data[:,b_index]==b,Eo_index] for b in range(bmin,bmax+1)]
          bands_cols.append(bands)
        elif type_calc == 'gw':
        #elif 'gw' in evalues:
          bands = [data[data[:,b_index]==b,E_index] for b in range(bmin,bmax+1)]
          bands_cols.append(bands)
        #elif 'corr' in evalues:
        #bands = [data[data[:,b_index]==b,DE_index] for b in range(bmin,bmax+1)]
        #bands_cols.append(bands)

        return kpoint_index, bands_cols

    def plot_qp_correction(self,tags=('qp',),lda=2,qp=3):
        if isinstance(tags,str): tags = (tags,)

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
        plt.plot()
        plt.show()

    def plot_gw(self,tags=('qp',),type_calc=('lda',),set_calc=(),rows=None):
        """
        Use this function to plot the quasiparticle energies from a GW calculation

        Arguments:
            type_calc: LDA, GW, or correction
            set_calc:  Plot a specific set, like EXXRLvcs. Default empty (all data)

        Example 1:
            a.plot_gw('qp', ('lda','gw'))

        Example 2:
            a.plot_gw('qp', ('lda','gw'),('EXXRLvcs',))
        """
        # The function is a bit redundant but it works
        '''
        Assign type of line
        '''
        line_dict  = {'lda':'-','gw':'--','corr':'-.'}
        '''
        Assign color with the number of calculations
        '''
        n_calculations = len(self.jsonfiles.keys())
        cmap = plt.get_cmap(self._colormap) #get color map
        colors = [cmap(i) for i in np.linspace(0, 1, n_calculations)]


        plot = False
        ax = plt.subplot(111)

        json_label = [ word[:-5] for word in self.jsonfiles.keys() ]


        # If we select a set of data
        if set_calc:
          for set_plot in set_calc:
            for i_json,json_filename in enumerate(sorted(self.jsonfiles.keys())):
              if set_plot in json_filename:
                for i_type,e_data in enumerate(type_calc):
                  kpoint_index, bands_cols = self.get_gw_bands(json_filename,e_data,json_filename[:-5])
                  label_data = '%s-%s' % (json_filename[:-5],e_data)
                  for i_bands,bands in enumerate(bands_cols):
                    for i_band,band in enumerate(bands): 
                      if i_band == 0:
                        ax.plot(kpoint_index,band,color=colors[i_json],linestyle=line_dict[e_data],label=label_data)
                      else:
                        ax.plot(kpoint_index,band,color=colors[i_json],linestyle=line_dict[e_data])
 
                        plot = True

        # If we plot all data
        if not set_calc: 
          for i_json,json_filename in enumerate(sorted(self.jsonfiles.keys())):
             for i_type,e_data in enumerate(type_calc):
               kpoint_index, bands_cols = self.get_gw_bands(json_filename,e_data,json_filename[:-5])
               label_data = '%s-%s' % (json_filename[:-5],e_data)

               for i_bands,bands in enumerate(bands_cols):
                 for i_band,band in enumerate(bands): 
                   if i_band == 0:
                     ax.plot(kpoint_index,band,color=colors[i_json],linestyle=line_dict[e_data],label=label_data)
                   else:
                     ax.plot(kpoint_index,band,color=colors[i_json],linestyle=line_dict[e_data])
 
                     plot = True

        if plot:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

            plt.title('GW quasiparticles on a mesh')
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':8})
            plt.show()

    def plot_bse(self,tags,cols=(2,),ax=None):
        """
        Use this function to plot the absorption spectrum calculated using the BSE
        cols: a list of indexes to select which columns from the file to plot

        Example:
            a.plot_bse('eps',cols=(2,))

            Will plot only files with 'eps' in their filename (absorption spectra)
            Will plot the second column (absorption spectra)
        """
        if ax is None:
            standalone = True
            ax = plt.gca()
        else:
            standalone = False
        plot = False

        colors = self.get_colors(tags)

        n=0
        for k in sorted(self.jsonfiles.keys()):
            for filename in list(self.jsonfiles[k]["data"].keys()):
                if all(i in filename for i in tags):
                    data = np.array( self.jsonfiles[k]["data"][filename] )

                    #select the color to plot with
                    color = colors[n]
                    n+=1

                    for col in cols:
                        x = data[:,0]
                        y = data[:,col-1]
                        label = filename.split('/')[-1]+" col=%d"%col
                        ax.plot(x,y,label=label,color=color)
                        plot = True
        if plot:
            ax.set_ylabel('Im$[\\chi(\omega)]$')
            ax.set_xlabel('$\omega$ (eV)')

            ax.legend(frameon=False)
            if standalone: plt.show()
        return ax

    def plot_spectral_function(self,tags):
        """
        Plot the spectral function
        """
        if isinstance(tags,str):
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
        for k in list(self.jsonfiles.keys()):
            if all(i in k for i in tags):
                print("\n%s"%k)
                for key,val in list(self.jsonfiles[k]["runtime"].items()):
                    print("%40s %10s %10s %10s"%(key,val[0],val[1],val[2]))

    def print_inputfiles(self):
        """
        Print all the inputfiles from all the json files
        """
        #iterate over the json files
        for k in list(self.jsonfiles.keys()):
            print("jsonfile: ", k)

            #get the jsonfile
            jsonfile = self.jsonfiles[k]

            for inputfile,content in list(jsonfile['inputfile'].items()):
                print("inputfile:", inputfile)
                y = YamboIn(filename=None)
                y.arguments = content["arguments"]
                y.variables = content["variables"]
                print(y)
                print()

    def __str__(self):
        s = ""
        for json_file in self.jsonfiles:
            s+="%s\n"%json_file
            for f in self.jsonfiles[json_file]['data']:
                s+="\t%s\n"%f
        return s
