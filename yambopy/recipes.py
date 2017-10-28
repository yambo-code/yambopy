#
# Copyright (C) 2017 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez, Alexandre Morlet, Fulvio Paleari
#
# All rights reserved.
#
# This file is part of yambopy
#
#
from yambopy import *
import os
from operator import itemgetter

#
# by Henrique Miranda. 
#
def pack_files_in_folder(folder,save_folder=None,mask='',verbose=True):
    """
    Pack the output files in a folder to json files
    """
    if not save_folder: save_folder = folder
    #pack the files in .json files
    for dirpath,dirnames,filenames in os.walk(folder):
        #check if the folder fits the mask
        if mask in dirpath:
            #check if there are some output files in the folder
            if ([ f for f in filenames if 'o-' in f ]):
                if verbose: print(dirpath)
                y = YamboOut(dirpath,save_folder=save_folder)
                y.pack()

#
# by Alejandro Molina-Sanchez
#
def breaking_symmetries(efield1,efield2=[0,0,0],folder='.',RmTimeRev=True):
    """
    Breaks the symmetries for a given field.
    Second field used in circular polarized pump configuration
    RmTimeRev : Remove time symmetry is set True by default
    """
    os.system('mkdir -p %s'%folder)
    os.system('cp -r database/SAVE %s'%folder)
    os.system('cd %s; yambo'%folder)
    ypp = YamboIn('ypp_ph -y -V all',folder=folder,filename='ypp.in')
    ypp['Efield1'] = efield1 # Field in the X-direction
    ypp['Efield2'] = efield2 # Field in the X-direction
    if RmTimeRev:
        ypp.arguments.append('RmTimeRev')   # Remove Time Symmetry
    ypp.write('%s/ypp.in'%folder)
    os.system('cd %s ; ypp_ph -F ypp.in'%folder )
    os.system('cd %s ; cd FixSymm; yambo '%folder )
    os.system('rm -r %s/SAVE'%folder)
    os.system('mv %s/FixSymm/SAVE %s/'%(folder,folder))
    os.system('rm -r %s/FixSymm'%folder)


#
# by Alexandre Morlet & Henrique Miranda
#

def analyse_gw(folder,var,bandc,kpointc,bandv,kpointv,pack,text,draw,verbose=False):
    """
    Study the convergence of GW calculations by looking at the change in band-gap value.

    The script reads from <folder> all results from <variable> calculations and display them.

    Use the band and k-point options (or change default values) according to the size of your k-grid and
    the location of the band extrema.
    """

    print('                  K-point   Band')
    print('Conduction state   %6d %6d'%(kpointc, bandc))
    print('   Valence state   %6d %6d'%(kpointv, bandv))

    # Packing results (o-* files) from the calculations into yambopy-friendly .json files
    if pack:
        if verbose: print('\nPacking...')
        pack_files_in_folder(folder,mask=var,verbose=verbose)
        pack_files_in_folder(folder,mask='reference',verbose=verbose)
        
    # importing data from .json files in <folder>
    if verbose: print('\nImporting data...')
    ya = YamboAnalyser(folder)
    
    # extract data according to relevant variable
    outvars = ya.get_data(tags=(var,'reference'))
    invars = ya.get_inputfiles_tag(var)
    tags = ya.get_tags(tags=(var,'reference'))

    # Get only files related to the convergence study of the variable,
    # ordered to have a smooth plot
    keys=[key for key in invars.keys() if key.startswith(var) or key=='reference.json']

    if len(keys) == 0: 
        raise ValueError('No files with this variable were found')

    if verbose:
        print('\nFiles detected:')
        for key in keys:
            print(" "*4,key)

    if verbose: print('\nComputing values...')
    ### Output

    # Unit of the variable :
    unit = invars[keys[0]]['variables'][var][1]

    # The following variables are used to make the script compatible with both short and extended output
    #kpindex = tags[keys[0]].tolist().index('K-point')
    kpindex = tags[keys[0]].tolist().index('Kpoint_index') # netcdf from json
    bdindex = tags[keys[0]].tolist().index('Band')
    e0index = tags[keys[0]].tolist().index('Eo')
    gwindex = tags[keys[0]].tolist().index('E-Eo')
    array = np.zeros((len(keys),2))

    for i,key in enumerate(keys):

        # input value
        # GbndRnge and BndsRnX_ are special cases
        if var.startswith('GbndRng') or var.startswith('BndsRnX'):
            # format : [1, nband, ...]
            array[i][0] = invars[key]['variables'][var][0][1]
        else:
            array[i][0] = invars[key]['variables'][var][0]

        # Output value (gap energy)
        # First the relevant lines are identified
        valence=[]
        conduction=[]
        for j in range(len(outvars[key]+1)):
            this_state = outvars[key][j]
            current_kpoint = this_state[kpindex]
            current_band   = this_state[bdindex]
            if   current_kpoint == kpointc and current_band == bandc:
                conduction = this_state
            elif current_kpoint == kpointv and current_band == bandv:
                valence = this_state
        # Then the gap can be calculated
        array[i][1] = conduction[e0index]+conduction[gwindex]-(valence[e0index]+valence[gwindex])

    #ascending order
    array = sorted( array, key=lambda x: x[0] )

    if text:
        os.system('mkdir -p analyse_%s'%folder)
        outname = './analyse_%s/%s_%s.dat'%(folder,folder,var)
        header = var+' ('+str(unit)+'), gap'
        np.savetxt(outname,array,delimiter='\t',header=header)
        if verbose: print('\nData saved to ',outname)

    if draw:
        plt.plot(array[:,0],array[:,1],'o-')
        plt.xlabel(var+' ('+unit+')')
        plt.ylabel('E_gw = E_lda + \Delta E')
        plt.savefig('%s.png'%var)
        if 'DISPLAY' in os.environ:
            plt.show()
    
    if verbose: print("\nDone")
#
# by Alexandre Morlet
#
def analyse_bse(folder,var,numbexc,intexc,degenexc,maxexc,pack,text,draw,verbose=False):
    """
    Using ypp, you can study the convergence of BSE calculations in 2 ways:
      Create a .png of all absorption spectra relevant to the variable you study
      Look at the eigenvalues of the first n "bright" excitons (given a threshold intensity)

    The script reads from <folder> all results from <variable> calculations for processing.
    The resulting pictures and data files are saved in the ./analyse_<folder>/ folder.

    Arguments:
        folder   -> Folder containing SAVE and convergence runs.
        var      -> Variable tested (e.g. FFTGvecs)
        numbexc  -> Number of excitons to read beyond threshold (default=2)
        intexc   -> Minimum intensity for excitons to be considered bright (default=0.05)
        degenexc -> Energy threshold under which different peaks are merged (eV) (default=0.01)
        maxexc   -> Energy threshold after which excitons are not read anymore (eV) (default=8.0)
        pack     -> Skips packing o- files into .json files (default: True)
        text     -> Skips writing the .dat file (default: True)
        draw     -> Skips drawing (plotting) the abs spectra (default: True)

    Returns:
        excitons -> energies of the first few excitons as funciton of some variable
        spectras -> absorption spectra for each variable

    """

    # Packing results (o-* files) from the calculations into yambopy-friendly .json files
    if pack: # True by default, False if -np used
        print('Packing ...')
        pack_files_in_folder(folder,mask=var)
        pack_files_in_folder(folder,mask='reference')

    # importing data from .json files in <folder>
    print('Importing data...')
    data = YamboAnalyser(folder)

    # extract data according to relevant var
    invars = data.get_inputfiles_tag(var)

    # Get only files related to the convergence study of the variable,
    # ordered to have a smooth plot
    keys=[key for key in invars.keys() if key.startswith(var) or key=='reference.json']

    if len(keys) == 0: 
        raise ValueError('No files with this variable were found')

    if verbose:
        print('Files detected:')
        for key in keys:
            print(key)

    # unit of the input value
    unit = invars[keys[0]]['variables'][var][1]

    ######################
    # Output-file filename
    ######################
    os.system('mkdir -p analyse_%s'%folder)
    outname = './analyse_%s/%s_%s'%(folder,folder,var)

    # Arrays that will contain the output
    excitons = []
    spectras = []

    # Loop over all calculations
    for key in keys:
        jobname=key.replace('.json','')
        print(jobname)

        # input value
        v = invars[key]['variables'][var][0]
        if type(v) == list:
            inp = v[1]
        else:
            inp = v

        print('Preparing JSON file. Calling ypp if necessary.')
        ### Creating the 'absorptionspectra.json' file
        # It will contain the exciton energies
        y = YamboOut(folder=folder,save_folder=folder)
        # Args : name of job, SAVE folder path, folder where job was run path
        a = YamboBSEAbsorptionSpectra(jobname,path=folder)
        # Get excitons values (runs ypp once)
        a.get_excitons(min_intensity=intexc,max_energy=maxexc,Degen_Step=degenexc)
        # Write .json file with spectra and eigenenergies
        a.write_json(filename=outname)

        ### Loading data from .json file
        with open(outname+'.json') as f:
            data = json.load(f)

        ### Plotting the absorption spectra
        spectras.append({'x': data['E/ev[1]'],
                         'y': data['EPS-Im[2]'],
                         'label': jobname})

        ### BSE spectra
        ### Axes : lines for exciton energies (disabled, would make a mess)
        #for n,exciton in enumerate(data['excitons']):
        #    plt.axvline(exciton['energy'])

        ### Creating array with exciton values (according to settings)
        l = [inp]
        for n,exciton in enumerate(data['excitons']):
            if n <= numbexc-1:
                l.append(exciton['energy'])
        excitons.append(l)
    
    #order in ascending order
    excitons,spectras = zip(*sorted(zip(excitons,spectras),key=lambda x: x[0][1]))

    if text:
        header = 'Columns : '+var+' (in '+unit+') and "bright" excitons eigenenergies in order.'

        ## Excitons energies
        #output on the screen
        print(header)
        for exc in excitons:
            x = exc[0]
            e = exc[1:]
            print("%8.4lf "%x+("%8.4lf"*len(e))%tuple(e))

        #save file
        filename = outname+'_excitons.dat'
        np.savetxt(filename,excitons,header=header)
        print(filename)

        ## Spectra
        filename = outname+'_spectra.dat'
        with open(filename,'w') as f:
            for spectra in spectras:
                label = spectra['label']
                f.write('#%s\n'%label)
                for x,y in zip(spectra['x'],spectra['y']):
                    f.write("%12.8e %12.8e\n"%(x,y))
                f.write('\n\n')
        print(filename)

    else:
        print('-nt flag : no text produced.')

    if draw:
        ## Exciton energy plots
        filename = outname+'_excitons.png'
        excitons = np.array(excitons)
        labels = [spectra['label'] for spectra in spectras]
        fig = plt.figure(figsize=(6,5))
        matplotlib.rcParams.update({'font.size': 15})
        plt.ylabel('1st exciton energy (eV)')
        plt.xticks(excitons[:,0],labels)
        plt.plot(excitons[:,0],excitons[:,1])
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        if 'DISPLAY' in os.environ:
            plt.show()
        print(filename)

        ## Spectra plots
        filename = outname+'_spectra.png'
        fig = plt.figure(figsize=(6,5))
        matplotlib.rcParams.update({'font.size': 15})
        for spectra in spectras:
            plt.plot(spectra['x'],spectra['y'],label=spectra['label'])
        plt.xlabel('$\omega$ (eV)')
        plt.ylabel('Im($\epsilon_M$)')
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        if 'DISPLAY' in os.environ:
            plt.show()
        print(filename)
    else:
        print('-nd flag : no plot produced.')

    print('Done.')

    return excitons, spectras

#
# by Fulvio Paleari & Henrique Miranda
#
def merge_qp(output,files,verbose=False):
    """
    Merge the quasiparticle databases produced by yambo
    """
    #read all the files and display main info in each of them
    print("=========input=========")
    filenames = [ f.name for f in files]
    datasets  = [ Dataset(filename) for filename in filenames]
    QP_table, QP_kpts, QP_E_E0_Z = [], [], []
    for d,filename in zip(datasets,filenames):
        _, nkpoints, nqps, _, nstrings = list(map(int,d['PARS'][:]))
        print("filename:    ", filename)
        if verbose:
            print("description:")
            for i in range(1,nstrings+1):
                print(''.join(d['DESC_strings_%05d'%i][0]))
        else:
            print("description:", ''.join(d['DESC_strings_%05d'%(nstrings)][0]))
        print()
        QP_table.append( d['QP_table'][:].T )
        QP_kpts.append( d['QP_kpts'][:].T )
        QP_E_E0_Z.append( d['QP_E_Eo_Z'][:] )

    # create the QP_table
    QP_table_save = np.vstack(QP_table)

    # create the kpoints table
    #create a list with the bigger size of QP_table
    nkpoints = int(max(QP_table_save[:,2]))
    QP_kpts_save = np.zeros([nkpoints,3])
    #iterate over the QP's and store the corresponding kpoint
    for qp_file,kpts in zip(QP_table,QP_kpts):
        #iterate over the kpoints and save the coordinates on the list
        for qp in qp_file:
            n1,n2,nk = list(map(int,qp))
            QP_kpts_save[nk-1] = kpts[nk-1]

    # create the QPs energies table
    QP_E_E0_Z_save = np.concatenate(QP_E_E0_Z,axis=1)

    #create reference file from one of the files
    netcdf_format = datasets[0].data_model
    fin  = datasets[0]
    fout = Dataset(output,'w',format=netcdf_format)

    variables_update = ['QP_table', 'QP_kpts', 'QP_E_Eo_Z']
    variables_save   = [QP_table_save.T, QP_kpts_save.T, QP_E_E0_Z_save]
    variables_dict   = dict(list(zip(variables_update,variables_save)))
    PARS_save = fin['PARS'][:]
    PARS_save[1:3] = nkpoints,len(QP_table_save)

    #create the description string
    kmin,kmax = np.amin(QP_table_save[:,2]),np.amax(QP_table_save[:,2])
    bmin,bmax = np.amin(QP_table_save[:,1]),np.amax(QP_table_save[:,1])
    description = "QP @ K %03d - %03d : b %03d - %03d"%(kmin,kmax,bmin,bmax)
    description_save = np.array([i for i in " %s"%description])

    #output data
    print("========output=========")
    print("filename:    ", output)
    print("description: ", description)

    #copy dimensions
    for dname, the_dim in list(fin.dimensions.items()):
        fout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

    #get dimensions
    def dimensions(array):
        return tuple([ 'D_%010d'%d for d in array.shape ])

    #create missing dimensions
    for v in variables_save:
        for dname,d in zip( dimensions(v),v.shape ):
            if dname not in list(fout.dimensions.keys()):
                fout.createDimension(dname, d)

    #copy variables
    for v_name, varin in list(fin.variables.items()):
        if v_name in variables_update:
            #get the variable
            merged = variables_dict[v_name]
            # create the variable
            outVar = fout.createVariable(v_name, varin.datatype, dimensions(merged))
            # Copy variable attributes
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            #save outvar
            outVar[:] = merged

        else:
            # create the variable
            outVar = fout.createVariable(v_name, varin.datatype, varin.dimensions)
            # Copy variable attributes
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            if v_name=='PARS':
                outVar[:] = PARS_save[:]
            elif v_name=='DESC_strings_%05d'%(nstrings):
                outVar[:] = varin[:]
                outVar[:,:len(description_save)] = description_save.T
            else:
                outVar[:] = varin[:]

    fout.close()


#
# by Alexandre Morlet, Fulvio Paleari & Henrique Miranda
#
def add_qp(output,add=[],substract=[],addimg=[],verbose=False):
    """
    Add quasiparticle lifetimes from multiple files
    """
    # Define filenames
    addf=[f.name for f in add]
    subf=[f.name for f in substract]
    addimgf=[f.name for f in addimg]
    filenames = addf+subf+addimgf

    if len(filenames) is 0:
        raise ValueError('No files passed to function.')


    # Init empty lists and dics
    sizes=[] # contains the various 'PARS'
    QP_table, QP_kpts, QP_E_E0_Z = {},{},{} # read value for each file
    qpdic = {} # used to calculate the final E (real part)
    qpdici = {} # used to calculate the final E (img part)

    # Read the files
    datasets  = [ Dataset(filename) for filename in filenames]

    print("\n    Reading input files\n")
    for d,f in zip(datasets,filenames):
        print("filename: %s"%f)
        # read sizes
        _, nkpoints, nqps, _, nstrings = list(map(int,d['PARS'][:]))
        sizes.append((f,(nkpoints,nqps,nstrings)))

        # Check if the number of kpoints is consistent
        # (Don't forget to break symmetries on every file for RT)
        if nkpoints!=sizes[0][1][0]:
            raise ValueError('File %s does not have the same number of kpoints'%f)

        # printing the description string
        # (breaking the symmetries doesn't update the descr)
        if verbose:
            print("description:")
            for i in range(1,nstrings+1):
                print(''.join(d['DESC_strings_%05d'%i][0]))
        else:
            print("description:", ''.join(d['DESC_strings_%05d'%(nstrings)][0]))

        # fill dictionaries with data for all files
        QP_table[f] = d['QP_table'][:].T
        QP_kpts[f]  = d['QP_kpts'][:].T
        QP_E_E0_Z[f]= d['QP_E_Eo_Z'][:]

        # Init qpdic & qpdici (going through each file in case the number of bands is different)
        # For qpdici, we assume Im(Eo)=0
        for (n1,n2,k),(E,Eo,Z) in zip(QP_table[f],QP_E_E0_Z[f][0]):
            qpdic[(n1,n2,k)]=Eo
            qpdici[(n1,n2,k)]=0

    print("Number of k points: %s\n"%nkpoints)

    # keys are sorted in the order yambo usually writes DBs
    qpkeys = sorted(list(qpdic.keys()),key=itemgetter(2,1))

    # For E, [0,:,:] is real part and [1,:,:] is img part
    QP_E_E0_Z_save = np.zeros((2,len(qpkeys),3))
    QP_table_save  = np.zeros((len(qpkeys),3))

    # create and init the QPs energies table

    # The E0 is simply written in the real part (is 0 in the img part) 
    # and Z = 1 (since we merge different calculation types)
    for i,(n1,n2,k) in enumerate(qpkeys):
        QP_E_E0_Z_save[0,i,1] = qpdic[(n1,n2,k)]
    QP_E_E0_Z_save[0,:,2] = 1
    QP_E_E0_Z_save[1,:,1] = 0
    QP_E_E0_Z_save[1,:,2] = 1


    # Add corrections in real part (-a files)
    for f in addf:
        print('Add E corr for real part :  %s'%f)
        for (n1,n2,k),(E,Eo,Z) in zip(QP_table[f],QP_E_E0_Z[f][0]):
            qpdic[(n1,n2,k)]+=E-Eo

    # Sub corrections in real part (-s files)
    for f in subf:
        print('Sub E corr for real part :  %s'%f)
        for (n1,n2,k),(E,Eo,Z) in zip(QP_table[f],QP_E_E0_Z[f][0]):
            qpdic[(n1,n2,k)]-=E-Eo

    # Add corrections in img part (-ai files)
    for f in addimgf:
        print('Add E corr for img part :  %s'%f)
        for (n1,n2,k),(E,Eo,Z) in zip(QP_table[f],QP_E_E0_Z[f][1]):
            qpdici[(n1,n2,k)]+=E-Eo


    # create the kpoints table
    # We put the restriction to have the same number of k points (same grid), so any file fits
    QP_kpts_save = QP_kpts[filenames[0]]

    # Filling the E column
    for i,(n1,n2,k) in enumerate(qpkeys):
        QP_table_save[i]=[n1,n2,k]
        QP_E_E0_Z_save[0,i,0]+=qpdic[(n1,n2,k)]
        QP_E_E0_Z_save[1,i,0]+=qpdici[(n1,n2,k)]


    ## Output file

    #create reference file from one of the files
    netcdf_format = datasets[0].data_model
    fin  = datasets[0]
    fout = Dataset(output,'w',format=netcdf_format)

    variables_update = ['QP_table', 'QP_kpts', 'QP_E_Eo_Z']
    variables_save   = [QP_table_save.T, QP_kpts_save.T, QP_E_E0_Z_save]
    variables_dict   = dict(list(zip(variables_update,variables_save)))
    PARS_save = fin['PARS'][:]
    PARS_save[1:3] = sizes[0][1][0],len(QP_table_save)

    #create the description string
    kmin,kmax = np.amin(QP_table_save[:,2]),np.amax(QP_table_save[:,2])
    bmin,bmax = np.amin(QP_table_save[:,1]),np.amax(QP_table_save[:,1])
    description = "QP @ K %03d - %03d : b %03d - %03d"%(kmin,kmax,bmin,bmax)
    description_save = np.array([i for i in " %s"%description])

    #output data
    print("\n    Producing output file\n")
    print("filename:    ", output)
    print("description: ", description)

    #copy dimensions
    for dname, the_dim in list(fin.dimensions.items()):
        fout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

    #get dimensions
    def dimensions(array):
        return tuple([ 'D_%010d'%d for d in array.shape ])

    #create missing dimensions
    for v in variables_save:
        for dname,d in zip( dimensions(v),v.shape ):
            if dname not in list(fout.dimensions.keys()):
                fout.createDimension(dname, d)

    #copy variables
    for v_name, varin in list(fin.variables.items()):
        if v_name in variables_update:
            #get the variable
            merged = variables_dict[v_name]
            # create the variable
            outVar = fout.createVariable(v_name, varin.datatype, dimensions(merged))
            # Copy variable attributes
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            #save outvar
            outVar[:] = merged

        else:
            # create the variable
            outVar = fout.createVariable(v_name, varin.datatype, varin.dimensions)
            # Copy variable attributes
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
            if v_name=='PARS':
                outVar[:] = PARS_save[:]
            elif v_name=='DESC_strings_%05d'%(nstrings):
                outVar[:] = varin[:]
                outVar[:,:len(description_save)] = description_save.T
            else:
                outVar[:] = varin[:]

    fout.close()



#
# by Henrique Miranda
#
def plot_excitons(filename,cut=0.2,size=20):
    from math import ceil, sqrt

    def get_var(dictionary,variables):
        """
        To have compatibility with different versions of yambo
        We provide a list of different possible tags
        """
        for var in variables:
            if var in dictionary:
                return dictionary[var]
        raise ValueError( 'Could not find the variables %s in the output file'%str(variables) )
    #
    # read file
    #
    f = open(filename)
    data = json.load(f)
    f.close()

    #
    # plot the absorption spectra
    #
    nexcitons = len(data['excitons'])
    print("nexitons", nexcitons)
    plt.plot(get_var(data,['E/ev','E/ev[1]']), get_var(data,['EPS-Im[2]' ]),label='BSE',lw=2)
    plt.plot(get_var(data,['E/ev','E/ev[1]']), get_var(data,['EPSo-Im[4]']),label='IP',lw=2)
    for n,exciton in enumerate(data['excitons']):
        plt.axvline(exciton['energy'])
    plt.xlabel('$\\omega$ (eV)')
    plt.ylabel('Intensity arb. units')
    plt.legend(frameon=False)
    plt.draw()

    #
    # plot excitons
    #

    #dimensions
    nx = int(ceil(sqrt(nexcitons)))
    ny = int(ceil(nexcitons*1.0/nx))
    print("cols:",nx)
    print("rows:",ny)
    cmap = plt.get_cmap("gist_heat_r")

    fig = plt.figure(figsize=(nx*3,ny*3))

    sorted_excitons = sorted(data['excitons'],key=lambda x: x['energy'])

    for n,exciton in enumerate(sorted_excitons):
        #get data
        w   = np.array(exciton['weights'])
        qpt = np.array(exciton['qpts'])

        #plot
        ax = plt.subplot(ny,nx,n+1)
        ax.scatter(qpt[:,0], qpt[:,1], s=size, c=w, marker='H', cmap=cmap, lw=0, label="%5.2lf (eV)"%exciton['energy'])
        ax.text(-cut*.9,-cut*.9,"%5.2lf (eV)"%exciton['energy'])

        # axis
        plt.xlim([-cut,cut])
        plt.ylim([-cut,cut])
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.set_aspect('equal')

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=0.01)

    #remove extension from file
    figure_filename = os.path.splitext(filename)[0]
    plt.savefig('%s.png'%figure_filename)
    if 'DISPLAY' in os.environ:
        plt.show()




