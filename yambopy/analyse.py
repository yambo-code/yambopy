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
from yambopy.plot import YambopyBandStructure
from yambopy.duck import isstring
from yambopy.lattice import red_car, rec_lat, expand_kpts, isbetween
from yambopy.io.inputfile import YamboIn

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

    def get_files_type(self,type,tags=None):
        """
        In all the json files present find the ones of a certain type
        the possible types are:
            report
            netcdf_gw
            netcdf_hf
        """
        json_files = self.jsonfiles
        files = {}

        #iterate over all the json files
        for json_filename in json_files.keys():
            json_file = json_files[json_filename]
            #all the output files in each json file
            for output_filename in json_file["files"]:
                output_file = json_file["files"][output_filename]
                if output_file["type"] == type:
                    filename, extension = os.path.splitext(json_filename)
                    files[filename] =  output_file

        #filter files with tags
        if tags:
            if isstring(tags):
                tags = (tags,)
            
            filenames = list(files.keys())
            for filename in filenames:
                if all(tag not in filename for tag in tags):
                    files.pop(filename)
        return files

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

    def get_gw_bands(self,tags=None,bs=None,type_calc=('ks','gw')):
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
        #get files that have a gw calculation in them
        gw_files = self.get_files_type('netcdf_gw',tags)

        #check if the dimensions of the files are consistent
        #TODO

        #create bandstructure class
        if not bs:
            bs = YambopyBandStructure()
 
        # add bandstructures of all the files
        for filename, content in gw_files.items():
            e0,e0imag = content['Eo']
            e,linewidths = content['E']
            ec,linewidths = content['E-Eo']
           
            #TODO move this section to YamboFileGW class
            #begin section            

            #get dimensions 
            band_index = np.array(content['Band'],dtype=int)
            band_min, band_max = min(band_index), max(band_index)
            nbands = band_max-band_min+1
            
            kpoint_index = np.array(content['Kpoint_index'],dtype=int)
            kpoint_min, kpoint_max = min(kpoint_index), max(kpoint_index)
            nkpoints = kpoint_max-kpoint_min+1

            #get arrys of bands and kpoints
            bands_e0 = np.zeros([nkpoints,nbands])
            bands_e  = np.zeros([nkpoints,nbands])
            for ei,e0i,ki,ni in zip(e,e0,kpoint_index,band_index):
                bands_e0[ki-1,ni-1] = e0i
                bands_e[ki-1,ni-1] = ei
            #end section

            #add bands
            if 'ks' in type_calc: 
                bs.add_bands(bands_e0,label=filename+' KS')
            if 'gw' in type_calc:
                bs.add_bands(bands_e, label=filename+' GW')

        return bs

    def get_gw_bands_path(self,path_label,tags=('qp',),type_calc=('lda',),rows=None):
        """ 
        Get the bands a path of k-points and find the points 
        in the regular mesh that correspond to points in the path
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

    def plot_ks(self,tags=None):
        """
        Use this function to plot the kohn sham energies from a GW calculation
        """
        #get bands from these files
        ks_bands = self.get_gw_bands(tags=tags,type_calc=('ks'))

        #plot the bands
        ks_bands.plot_show()
        return ks_bands

    def plot_gw(self,tags=None):
        """
        Use this function to plot the quasiparticle energies from a GW calculation
        """
        #get bands from these files
        gw_bands = self.get_gw_bands(tags=tags,type_calc=('gw'))

        #plot the bands
        gw_bands.plot_show()
        return gw_bands

    def plot_gw_path(self,path_label,tags=('qp',),type_calc=('lda',),set_calc=(),rows=None):
        """
        Create a path of k-points and find the points in the 
        regular mesh that correspond to points in the path
        Use these points to plot the GW band structure.
        """
        # assign type of line
        line_dict  = {'lda':'-','gw':'--','corr':'-.'}
        # assign color with the number of calculations
        n_calculations = len(self.jsonfiles.keys())
        cmap = plt.get_cmap(self._colormap) #get color map
        colors = [cmap(i) for i in np.linspace(0, 1, n_calculations)]

        path = np.array([p[0] for p in path_label])
        labels = [p[1] for p in path_label]
        plot = False

        fig = plt.figure()
        ax = plt.subplot(111)

        #implemente HERE!!

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
                print(y+'\n')

    def __str__(self):
        s = ""
        for json_file in self.jsonfiles:
            s+="%s\n"%json_file
            for f in self.jsonfiles[json_file]['data']:
                s+="\t%s\n"%f
        return s
