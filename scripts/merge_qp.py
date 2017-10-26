from __future__ import print_function
# Copyright (C) 2016 Fulvio Paleari, Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from builtins import zip
from builtins import map
from builtins import range
from yamboparser import *
from os import *
import argparse

def merge_qp(output,files,verbose=False):
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
    fin  = datasets[0]
    fout = Dataset(output,'w') 

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Join different NetCDF quasi-particle databases')
    parser.add_argument('files', nargs='+', type=argparse.FileType('r'))
    parser.add_argument('-o','--output',                       help='Output filename', default='ndb_out.QP')
    parser.add_argument('-v','--verbose', action="store_true", help='Verbose mode')
    args = parser.parse_args()

    if args.files is None:
        parser.print_help()
        exit()

    output  = args.output
    files   = args.files
    verbose = args.verbose
    merge_qp(output,files,verbose)




