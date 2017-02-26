# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez, Alexandre Morlet
# All rights reserved.
#
# This file is part of yambopy
#
#
from yambopy import *
import os

#
# by Henrique Miranda
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
                print dirpath
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
# by Alexandre Morlet
#
def analyse_bse( folder, var, numbexc, intexc, degenexc, maxexc, pack=True, text=True, draw=True ):
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
    """
            
    # Packing results (o-* files) from the calculations into yambopy-friendly .json files
    if pack: # True by default, False if -np used
        print 'Packing ...'
        pack_files_in_folder(folder,mask=var)
        pack_files_in_folder(folder,mask='reference')
        print 'Packing done.'
    else:
        print 'Packing skipped.'

    # importing data from .json files in <folder>
    print 'Importing...'
    data = YamboAnalyser(folder)

    # extract data according to relevant var
    invars = data.get_inputfiles_tag(var)

    # Get only files related to the convergence study of the variable,
    # ordered to have a smooth plot
    keys=[]
    sorted_invars = sorted(invars.items(), key=operator.itemgetter(1))

    for i in range(0,len(sorted_invars)):
        key=sorted_invars[i][0]
        if key.startswith(var) or key=='reference.json':
            keys.append(key)
    print 'Files detected: ',keys

    # unit of the input value
    unit = invars[keys[0]]['variables'][var][1]

    ######################
    # Output-file filename
    ######################
    os.system('mkdir -p analyse_%s'%folder)
    outname = './analyse_%s/%s_%s'%(folder,folder,var)

    # Array that will contain the output
    excitons = []

    # Loop over all calculations
    for key in keys:
        jobname=key.replace('.json','')
        print jobname

        # input value
        # BndsRn__ is a special case
        if var.startswith('BndsRnX'):
        # format : [1, nband, ...]
            inp = invars[key]['variables'][var][0][1]
        else:
            inp = invars[key]['variables'][var][0]

        print 'Preparing JSON file. Calling ypp if necessary.'
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
        f = open(outname+'.json')
        data = json.load(f)
        f.close()
        print 'JSON file prepared and loaded.'

        ### Plotting the absorption spectra
        # BSE spectra
        plt.plot(data['E/ev[1]'], data['EPS-Im[2]'],label=jobname,lw=2)
    #   # Axes : lines for exciton energies (disabled, would make a mess)
    #   for n,exciton in enumerate(data['excitons']):
    #       plt.axvline(exciton['energy'])

        ### Creating array with exciton values (according to settings)
        l = [inp]
        for n,exciton in enumerate(data['excitons']):
            if n <= numbexc-1:
                l.append(exciton['energy'])

        excitons.append(l)

    if text:
        header = 'Columns : '+var+' (in '+unit+') and "bright" excitons eigenenergies in order.'
        print excitons
        np.savetxt(outname+'.dat',excitons,header=header)
        #np.savetxt(outname,excitons,header=header,fmt='%1f')
        print outname+'.dat'
    else:
        print '-nt flag : no text produced.'

    if draw:
        plt.xlabel('$\omega$ (eV)')
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.legend()
        #plt.draw()
        #plt.show()
        plt.savefig(outname+'.png', bbox_inches='tight')
        print outname+'.png'
    else:
        print '-nd flag : no plot produced.'

    print 'Done.'

