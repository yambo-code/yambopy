# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function, division
import os
import json
import re
import numpy as np
import matplotlib.pyplot as plt
from yambopy.plot.bandstructure import YambopyBandStructure
from yambopy.kpoints import get_path
from yambopy.tools.duck import isstring
from yambopy.io.inputfile import YamboIn
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.plot.plotting import add_fig_kwargs
from yambopy.tools.string import marquee

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
        !!!!!!
        !!!!!! This function was buggy. Let's check it better
        !!!!!!
        Select the colors according to the number of files to plot
        the files to plot are the ones that have all the tags in their name
        """
        import matplotlib.pyplot as plt
        #count the number of files
        nfiles = 0
        for k in self.jsonfiles.keys():
            for filename in list(self.jsonfiles[k]['files'].keys()):
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

    def get_bands(self,tags=None,path_kpoints=None,type_calc=('ks','gw')):
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
        ks_bandstructure, qp_bandstructure = None, None

        # add bandstructures of all the files
        for filename, content in gw_files.items():
            e0,e0imag = content['Eo']
            e,linewidths = content['E']
            ec,linewidths = content['E-Eo']
            kpoints = content['Kpoint']

            #TODO move this section to YamboFileGW class
            #begin section            

            #get dimensions 
            band_index = np.array(content['Band'],dtype=int)
            band_min, band_max = min(band_index), max(band_index)
            nbands = band_max-band_min+1
           
            kpoint_index = np.array(content['Kpoint_index'],dtype=int)
            kpoint_min, kpoint_max = min(kpoint_index), max(kpoint_index)
            nkpoints = kpoint_max-kpoint_min+1

            #get arrays of bands and kpoints
            bands_e0 = np.zeros([nkpoints,nbands])
            bands_e  = np.zeros([nkpoints,nbands])
            for ei,e0i,ki,ni in zip(e,e0,kpoint_index,band_index):
                nkpoint = ki-kpoint_min
                nband = ni-band_min
                bands_e0[nkpoint,nband] = e0i
                bands_e[nkpoint,nband] = ei
            #end section

            if path_kpoints:
                #get data from json file
                jsonfile = list(self.jsonfiles.values())[0]
                lat = YamboLatticeDB.from_dict(jsonfile['lattice'])
                kpoints, bands_indexes, path_car = get_path(lat.car_kpoints,lat.rlat,None,path_kpoints)  
                bands_e0 = bands_e0[bands_indexes]
                bands_e  = bands_e[bands_indexes] 

            #add bands
            if 'ks' in type_calc: 
                ks_bandstructure = YambopyBandStructure(bands_e0,kpoints,label=filename)
            if 'gw' in type_calc:
                qp_bandstructure = YambopyBandStructure(bands_e, kpoints,label=filename)

        return ks_bandstructure, qp_bandstructure

    @add_fig_kwargs
    def plot_ks(self,path=None,tags=None):
        """
        Use this function to plot the kohn sham energies from a GW calculation
        """
        #get bands from these files
        ks_bands = self.get_bands(tags=tags,path=path,type_calc=('ks'))[0]
        
        #plot the bands
        return ks_bands.plot(show=False)

    @add_fig_kwargs
    def plot_gw(self,path_kpoints=None,tags=None,**kwargs):
        """
        Use this function to plot the quasiparticle energies from a GW calculation
        """
        print('tags')
        print(tags)
        print()
        #get bands from these files
        gw_bands = self.get_bands(tags=tags,path_kpoints=path_kpoints,type_calc=('gw',))[1]

        #plot the bands
        return gw_bands.plot(show=False)

    def plot_bse(self,tags,cols=(2,),ax=None,png_file=False):
        """
        Use this function to plot the absorption spectrum calculated using the BSE
        cols: a list of indexes to select which columns from the file to plot

        Example:
            a.plot_bse('eps_q1',cols=(2,))

            Will plot only files with 'eps_q1' in their filename (absorption spectra)
            Will plot the second column (absorption spectra)


            a.plot_bse(('eps_q1,'FFTGvecs'),cols=(2,))
            Will plot only files with 'eps_q1' in their filename (absorption spectra)
            and for FFTGvecs variable

        !!!!!
        !!!!! Problems here. I have removed there reference to data ["data"]

        """
        import matplotlib.pyplot as plt
        if ax is None:
            standalone = True
            ax = plt.gca()
        else:
            standalone = False
        plot = False

        colors = self.get_colors(tags)

        n=0
        for k in sorted(self.jsonfiles.keys()):
            for filename in list(self.jsonfiles[k]["files"].keys()):
                if all(i in filename for i in tags):
                    # I prefer to work directly with the dictionary...
                    #data = np.array( self.jsonfiles[k]["files"][filename] )
                    data = self.jsonfiles[k]["files"][filename]
                    #select the color to plot with
                    color = colors[n]
                    n+=1

                    for col in cols:
                        #x = data[:,0]
                        #y = data[:,col-1]
                        x = data['E/ev[1]']
                        y = data['EPS-Im[2]']
                        label = filename.split('/')[-1]+" col=%d"%col
                        # Should we clean the label?
                        label=label.replace('.eps_q1_haydock_bse','')
                        label=label.replace('o-','')
                        ax.plot(x,y,label=label,color=color)
                        plot = True
        if plot:
            ax.set_ylabel('Im$[\\chi(\omega)]$')
            ax.set_xlabel('$\omega$ (eV)')

            ax.legend(frameon=False,loc=1)
            if standalone: plt.show()
        if png_file:
            plt.savefig('%s.png' % label[:8])

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

    def plot_gw_all_kpoints_convergence(self,tag=None):
        '''
        Function to plot the GW-QPs energies of all k-points
        for a given tag
        Please test
        '''
        import matplotlib.pyplot as plt
        # 1. Find the json files with the tag given: FFTGvecs, BndsRnXp, etc.

        tag_list = []

        for word in self.jsonfiles.keys():
            if tag in word:
               tag_list.append(word.replace('.json','') )

        ntags = len(tag_list)

        cmap = plt.get_cmap('rainbow')
        colors=[cmap(i) for i in np.linspace(0,1,ntags)]

        bands_tag = []

        # 2. Get the bands of all keys
        for it in range(ntags):
            bands_tag.append(self.get_bands(tags=tag_list[it],type_calc=('gw',))[1])
            
        # 3. Get the bands of all keys
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for it in range(ntags):
            bands_tag[it].plot_ax(ax,color_bands=colors[it],c_label=tag_list[it],legend=True)

        plt.show()

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        for json_file in self.jsonfiles:
            app("%s"%json_file)
            for f in self.jsonfiles[json_file]['files']:
                app("\t%s"%f)
        return "\n".join(lines)
