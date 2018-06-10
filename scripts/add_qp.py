# Copyright (C) 2018 Alexandre Morlet, Fulvio Paleari, Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
from __future__ import print_function
from builtins import zip
from builtins import map
from builtins import range
from yamboparser import *
from os import *
import argparse
from operator import itemgetter

def add_qp(output,add=[],substract=[],addimg=[],verbose=False):
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
        print("filename:    ", f)
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

    # The E0 is simply written in the real part (is 0 in the img part) and Z = 1 (since we merge different calculation types)
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
    for dname, the_dim in fin.dimensions.items():
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
    for v_name, varin in fin.variables.items():
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Join different NetCDF quasi-particle databases')
    parser.add_argument('-a','--add', nargs='+', type=argparse.FileType('r'), help="Add the real part to the final db",default=[])
    parser.add_argument('-s','--substract', nargs='+', type=argparse.FileType('r'), help="Substract the real part to the final db", default=[])
    parser.add_argument('-ai','--addimg', nargs='+', type=argparse.FileType('r'), help="Add the imaginary part to the final db",default=[])
    parser.add_argument('-o','--output',                       help='Output filename', default='ndb_out.QP')
    parser.add_argument('-v','--verbose', action="store_true", help='Verbose mode')
    args = parser.parse_args()


    output  = args.output
    add     = args.add
    substract = args.substract
    addimg  = args.addimg
    verbose = args.verbose
    add_qp(output,add,substract,addimg,verbose)
