# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez, Alexandre Morlet, Fulvio Paleari
#
# All rights reserved.
#
# This file is part of yambopy
#
#
import os
from operator import itemgetter
from collections import OrderedDict
from yambopy import *

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
    ypp = YamboIn.from_runlevel('-y -V all',executable='ypp',folder=folder,filename='ypp.in')
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

    #find all ndb.QP files in the folder
    io = OrderedDict()
    for root, dirs, files in os.walk(folder):
        #get starting name of folder
        basename = os.path.basename(root)
        #look into folders starting with var or reference
        if any( [basename.startswith(v) for v in [var,'reference']] ):
            for filename in files:
                if filename != 'ndb.QP': continue
                #get ndb.QP file in folder
                io[basename] = ( YamboIn.from_file(folder=folder,filename="%s.in"%basename),
                                 YamboQPDB.from_db(folder=root,filename=filename) )

    #consistency check
    #TODO


    convergence_data = []

    for basename, (inp,out) in io.items():
        #get input
        value, unit = inp[var]

        #get qp value
        # Be careful because the array of eigenvalues is defined now in another way
        eigenvalues_dft, eigenvalues_qp, lifetimes, z = out.get_qps()
        #save result
        qp_gap = eigenvalues_qp[kpointc-out.min_kpoint,bandc-out.min_band] - eigenvalues_qp[kpointv-out.min_kpoint,bandv-out.min_band]

        #check type of variable
        if isinstance(value,list): value = value[1]
        convergence_data.append([value,qp_gap])

    convergence_data = np.array(sorted(convergence_data))
    if convergence_data.dtype == 'object': raise ValueError('Unknown type of variable')

    if text:
        output_folder = 'analyse_%s'%folder
        if not os.path.isdir(output_folder): os.mkdir(output_folder)
        outname = os.path.join(output_folder,'%s_%s.dat'%(folder,var))
        header = var+' ('+str(unit)+')'
        np.savetxt(outname,convergence_data,delimiter='\t',header=header)

    if draw:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(convergence_data[:,0],convergence_data[:,1],'o-')
        ax.set_xlabel(var+' ('+unit+')')
        ax.set_ylabel('$E_{GW} = E_{LDA} + \Delta E$')
        plt.show()
        fig.savefig('%s.png'%var)

#
# by Alexandre Morlet
#
def analyse_bse(folder,var,numbexc,intexc,degenexc,maxexc,text,draw,verbose=False):
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
        text     -> Skips writing the .dat file (default: True)
        draw     -> Skips drawing (plotting) the abs spectra (default: True)
    """

    #find the save folder
    lat = YamboSaveDB.from_db_file(os.path.join(folder,'SAVE'))

    #find all ndb.BS_diago_Q01 files in the folder
    io = OrderedDict()
    for root, dirs, files in os.walk(folder):
        #get starting name of folder
        basename = os.path.basename(root)
        #look into folders starting with var or reference
        if any( [basename.startswith(v) for v in [var,'reference']] ):
            for filename in files:
                if filename != 'ndb.BS_diago_Q01': continue
                #get ndb.BS_diago_Q01 file in folder
                io[basename] = ( YamboIn.from_file(folder=folder,filename="%s.in"%basename),
                                 YamboExcitonDB.from_db_file(lat,folder=root,filename=filename) )

    #TODO consistency check
    exciton_energies = []
    exciton_spectras = []
    for basename, (inp,out) in io.items():
        #get input
        value, unit = inp[var]

        #get exiton energies
        exciton_energy = out.eigenvalues.real

        #get excitonic spectra
        exciton_spectra = out.get_chi()

        #check type of variable
        if isinstance(value,list): value = value[1]
        exciton_energies.append([value,exciton_energy])
        exciton_spectras.append([value,exciton_spectra])

    exciton_spectras = sorted(exciton_spectras,key=lambda x: x[0])
    exciton_energies = sorted(exciton_energies,key=lambda x: x[0])

    #save a file with the exciton eneergies
    output_folder = 'analyse_%s'%folder
    if not os.path.isdir(output_folder): os.mkdir(output_folder)
    output_file = '%s_exciton_energies.dat'%var
    with open(os.path.join(output_folder,output_file),'w') as f:
        header = "%s (%s)\n"%(var,unit) if unit else "%s\n"%var
        f.write(header)
        for value,energies in exciton_energies:
            f.write("{} ".format(value)+("%10.6lf "*numbexc)%tuple(energies[:numbexc])+"\n")

    import matplotlib.pyplot as plt
    ## Exciton spectra plots
    filename = 'exciton_spectra.png'
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(1,1,1)

    #plot the spectra
    cmap = plt.get_cmap('viridis')
    nspectra = len(exciton_spectras)
    for i,(value,(w,spectra)) in enumerate(exciton_spectras):
        plt.plot(w,spectra.imag,c=cmap(i/nspectra),label="{} = {} {}".format(var,value,unit))

    ## Spectra plots
    ax.set_xlabel('$\omega$ (eV)')
    ax.set_ylabel('Im($\epsilon_M$)')
    ax.legend(frameon=False)
    output_file = '%s_exciton_spectra.pdf'%var
    fig.savefig(os.path.join(output_folder,output_file))
    if draw: plt.show()

#
# by Fulvio Paleari & Henrique Miranda
#
def merge_qp(output,files):
    """
    Merge the quasiparticle databases produced by yambo
    """
    #read all the files and display main info in each of them
    filenames = [ f.name for f in files]
    datasets  = [ Dataset(filename) for filename in filenames]
    #call compatibility version if old dataset detected 
    try: 
        qp_test = datasets[0]['QP_E']
    except IndexError:    
        print("Old version of database detected. Switching to compatibility mode.")
        try: qp_test = datasets[0]['QP_E_Eo_Z']
        except IndexError: raise IndexError('Problem with the databases')
        finally: merge_qp_compatibility(output,files,verbosity=False)
    finally:

        print("=========input=========")
        QP_table, QP_kpts, QP_E, QP_E0, QP_Z = [], [], [], [], []
        for d,filename in zip(datasets,filenames):
            pars_valid = [ par for par in d['PARS'][:] if not  np.ma.is_masked(par) ] # Fix to exclude empty elements in database list (masked by default by python)
            PARS = list(map(int,pars_valid))
            nkpoints, nqps, nstrings = PARS[1],PARS[2],PARS[-1]
            #_, nkpoints, nqps, _, nstrings = list(map(int,d['PARS'][:]))
            print("filename:    ", filename)
            QP_table.append( d['QP_table'][:].T )
            QP_kpts.append( d['QP_kpts'][:].T )
            QP_E.append( d['QP_E'][:] )
            QP_E0.append( d['QP_Eo'][:] )
            QP_Z.append( d['QP_Z'][:] )

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
                try:               n1,n2,nk = list(map(int,qp))
                except ValueError: n1,n2,nk,ns = list(map(int,qp))
                QP_kpts_save[nk-1] = kpts[nk-1]

        # create the QPs energies table
        QP_E_save  = np.concatenate(QP_E,axis=0)
        QP_E0_save = np.concatenate(QP_E0)
        QP_Z_save  = np.concatenate(QP_Z,axis=0)

        #create reference file from one of the files
        netcdf_format = datasets[0].data_model
        fin  = datasets[0]
        fout = Dataset(output,'w',format=netcdf_format)

        variables_update = ['QP_table', 'QP_kpts', 'QP_E', 'QP_Eo', 'QP_Z']
        variables_save   = [QP_table_save.T, QP_kpts_save.T, QP_E_save, QP_E0_save, QP_Z_save]
        variables_dict   = dict(list(zip(variables_update,variables_save)))
        PARS_save = fin['PARS'][:]
        PARS_save[1:3] = nkpoints,len(QP_table_save)

        #create the description string
        kmin,kmax = np.amin(QP_table_save[:,2]),np.amax(QP_table_save[:,2])
        bmin,bmax = np.amin(QP_table_save[:,1]),np.amax(QP_table_save[:,1])
        description = "QP @ K %03d - %03d : b %03d - %03d"%(kmin,kmax,bmin,bmax)
        description_save = np.array([i for i in " %s"%description])
        QP_k_range, QP_b_range = [kmin,kmax], [bmin,bmax]


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
                elif v_name=='QP_QP_@_state_1_K_range':
                    outVar[:]=QP_k_range
                elif v_name=='QP_QP_@_state_1_b_range':
                    outVar[:]=QP_b_range
                else:
                    outVar[:] = varin[:]

        fout.close()


#
# Old version by Alexandre Morlet, Fulvio Paleari & Henrique Miranda
# Updated by Daniel Murphy
#
def add_qp(output,add=[],subtract=[],addimg=[],verbose=False):
    """
    Add quasiparticle lifetimes from multiple files
    """
    # Define filenames
    addf=[f.name for f in add]
    subf=[f.name for f in subtract]
    addimgf=[f.name for f in addimg]
    filenames = addf+subf+addimgf
    # Read the files
    datasets  = [ Dataset(filename) for filename in filenames]

    if len(filenames) == 0:
        raise ValueError('No files passed to function.')

    #call compatibility version if old dataset detected
    try:
        qp_test = datasets[0]['QP_E']
    except IndexError:
        print("Old version of database detected. Switching to compatibility mode.")
        try: qp_test = datasets[0]['QP_E_Eo_Z']
        except IndexError: raise IndexError('Problem with the databases')
        finally: add_qp_compatibility(output,add,subtract,addimg,verbose=False)
    finally:

        # Init empty lists and dics
        sizes=[] # contains the various 'PARS'
        QP_table, QP_kpts, QP_E, QP_E0, QP_Z = {},{},{},{},{} # read value for each file
        qpdic = {} # used to calculate the final E (real part)
        qpdici = {} # used to calculate the final E (img part)

        print("\n    Reading input files\n")
        for d,f in zip(datasets,filenames):
            print("filename: %s"%f)
            # read sizes
            PARS = list(map(int, d["PARS"][:]))
            nkpoints, nqps, nstrings = PARS[1], PARS[2], PARS[-1]
            sizes.append((f,(nkpoints,nqps,nstrings)))

            # Check if the number of kpoints is consistent
            # (Don't forget to break symmetries on every file for RT)
            if nkpoints!=sizes[0][1][0]:
                raise ValueError('File %s does not have the same number of kpoints'%f)

            # So far, no verbose option
            #if verbose:

            # fill dictionaries with data for all files
            QP_table[f] = d['QP_table'][:].T
            QP_kpts[f]  = d['QP_kpts'][:].T
            QP_E[f] = d["QP_E"][:]
            QP_E0[f] = d["QP_Eo"][:]
            QP_Z[f] = d["QP_Z"][:]

            # Init qpdic & qpdici (going through each file in case the number of bands is different)
            # We assume Im(E0)=0
            for (n1,n2,k),(E0) in zip(QP_table[f],QP_E0[f]):
                qpdic[(n1,n2,k)]=E0
                qpdici[(n1,n2,k)]=0.

        print("Number of k points: %s\n"%nkpoints)

        # keys are sorted in the order yambo usually writes DBs
        # This is: [ (i_b1,i_b2=i_b1,i_k) ... ] looping first on i_b1
        qpkeys = sorted(list(qpdic.keys()),key=itemgetter(2,1))

        # For E, [:,0] is real part and [:,1] is img part
        QP_table_save = np.zeros((len(qpkeys), 3))
        QP_E_save = np.zeros((len(qpkeys), 2))
        QP_E0_save = np.zeros((len(qpkeys), 1))
        QP_Z_save = np.ones((len(qpkeys), 1))

        # Add corrections in real part (-a files)
        for f in addf:
            print('Add E corr for real part :  %s'%f)
            for (n1,n2,k),(E),(Eo) in zip(QP_table[f],QP_E[f][:,0],QP_E0[f][:]):
                qpdic[(n1,n2,k)]+=E-Eo

        # Sub corrections in real part (-s files)
        for f in subf:
            print('Sub E corr for real part :  %s'%f)
            for (n1,n2,k),(E),(Eo) in zip(QP_table[f],QP_E[f][:,0],QP_E0[f][:]):
                qpdic[(n1,n2,k)]-=E-Eo

        # Add corrections in img part (-ai files)
        for f in addimgf:
            print('Add E corr for img part :  %s'%f)
            for (n1,n2,k),(E) in zip(QP_table[f],QP_E[f][:,1]):
                qpdici[(n1,n2,k)]+=E

        # create the kpoints table
        # We put the restriction to have the same number of k points (same grid), so any file fits
        QP_kpts_save = QP_kpts[filenames[0]]

        # Filling the E column
        for i,(n1,n2,k) in enumerate(qpkeys):
            QP_table_save[i]=[n1,n2,k]
            QP_E_save[i, 0] += qpdic[(n1, n2, k)]
            QP_E_save[i, 1] += qpdici[(n1, n2, k)]
        QP_E0_save = array(QP_E0)

        ## Output file

        #create reference file from one of the files
        netcdf_format = datasets[0].data_model
        fin  = datasets[0]
        fout = Dataset(output,'w',format=netcdf_format)

        variables_update = ['QP_table','QP_kpts','QP_E','QP_E0','QP_Z']
        variables_save   = [QP_table_save.T, QP_kpts_save.T, QP_E_save,QP_E0_save,QP_Z_save]
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

def merge_qp_compatibility(output,files,verbose=False):
    """
    Merge the quasiparticle databases produced by yambo if they are in the old format
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
def add_qp_compatibility(output,add=[],subtract=[],addimg=[],verbose=False):
    """
    Add quasiparticle lifetimes from multiple files in the old format
    """
    # Define filenames
    addf=[f.name for f in add]
    subf=[f.name for f in subtract]
    addimgf=[f.name for f in addimg]
    filenames = addf+subf+addimgf

    if len(filenames) == 0:
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
