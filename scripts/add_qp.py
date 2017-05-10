# Copyright (C) 2017 Alexandre Morlet, Fulvio Paleari, Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from yamboparser import *
from os import *
import argparse

def add_qp(output,add,substract=[],addimg=[],verbose=False):
    # Define filenames
    addf=[f.name for f in add]
    subf=[f.name for f in substract]
    addimgf=[f.name for f in addimg]
    filenames = addf+subf+addimgf

    # Init empty lists and dics
    sizes=[] # contains the various 'PARS'
    QP_table, QP_kpts, QP_E_E0_Z = {},{},{} # read value for each file
    qpdic = {} # used to calculate the final E

    # Read the files
    datasets  = [ Dataset(filename) for filename in filenames]
    netcdf_format = datasets[0].data_model
    for d,filename in zip(datasets,filenames):
        # read sizes
        _, nkpoints, nqps, _, nstrings = map(int,d['PARS'][:])
        sizes.append((filename,(nkpoints,nqps,nstrings)))
        # Check if the number of kpoints is the same
        if nkpoints!=sizes[0][1][0]:
            print 'File %s does not have the same number of kpoints'%filename
            exit()
        print "filename:    ", filename

        # printing the description string
        # NB : breaking the symmetries don't update the descr
        if verbose:
            print "description:"
            for i in xrange(1,nstrings+1):
                print ''.join(d['DESC_strings_%05d'%i][0])
        else:
            print "description:", ''.join(d['DESC_strings_%05d'%(nstrings)][0])
        print "Number of k points: ", nkpoints

        # fill dictionaries with data for all files
        QP_table[filename] = d['QP_table'][:].T 
        QP_kpts[filename]  = d['QP_kpts'][:].T 
        QP_E_E0_Z[filename]= d['QP_E_Eo_Z'][:] 
        # Init qpdic (going through each file in case the number of bands is different)
        for (n1,n2,k) in QP_table[filename]:
            qpdic[(n1,n2,k)]=0
    
    # Add what needs to be added (-a files)
    for filename in addf:
        for (n1,n2,k),(E,Eo,Z) in zip(QP_table[filename],QP_E_E0_Z[filename][0]):
            qpdic[(n1,n2,k)]+=E

    # Sub what needs to be substracted (-s files)
    for filename in subf:
        for (n1,n2,k),(E,Eo,Z) in zip(QP_table[filename],QP_E_E0_Z[filename][0]):
            qpdic[(n1,n2,k)]-=E

    #print qpdic
    #exit()
#    for (f,s) in sizes:
#        print f
#        print s
#        if f in addf: print 'True'
#    exit()



    # create the QP_table
#    QP_table_save = np.vstack(QP_table)

    # create the kpoints table
    # We put the restriction to have the same number of k points (same grid), so any file fits
    QP_kpts_save = QP_kpts[filenames[0]]


    # create the QPs energies table
    qpkeys = sorted(qpdic.keys())

    QP_E_E0_Z_save = np.zeros((len(qpkeys),3))
    QP_table_save  = np.zeros((len(qpkeys),3))
    # The E0 and Z are taken from the first file, as E0 should be the same everywhere and we have no use for Z
    QP_E_E0_Z_save[:,1] = QP_E_E0_Z[filenames[0]][0][:,1]
    QP_E_E0_Z_save[:,2] = QP_E_E0_Z[filenames[0]][0][:,2]

    # Filling the E column
    for i,(n1,n2,k) in enumerate(qpkeys):
        print n1,n2,k
        QP_table_save[i]=[n1,n2,k]
        QP_E_E0_Z_save[i,0]=qpdic[(n1,n2,k)]
    
    print QP_table_save
    print QP_E_E0_Z_save
    # Ok 'til here


    #create reference file from one of the files
    fin  = datasets[0]
    fout = Dataset(output,'w',format=netcdf_format) 

    variables_update = ['QP_table', 'QP_kpts', 'QP_E_Eo_Z']
    variables_save   = [QP_table_save.T, QP_kpts_save.T, QP_E_E0_Z_save]
    variables_dict   = dict(zip(variables_update,variables_save)) 
    PARS_save = fin['PARS'][:]
    PARS_save[1:3] = nkpoints,len(QP_table_save)

    #create the description string
    kmin,kmax = np.amin(QP_table_save[:,2]),np.amax(QP_table_save[:,2])
    bmin,bmax = np.amin(QP_table_save[:,1]),np.amax(QP_table_save[:,1])
    description = "QP @ K %03d - %03d : b %03d - %03d"%(kmin,kmax,bmin,bmax)
    description_save = np.array([i for i in " %s"%description])

    #output data
    print "========output========="
    print "filename:    ", output
    print "description: ", description

    #copy dimensions
    for dname, the_dim in fin.dimensions.iteritems():
        fout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

    #get dimensions
    def dimensions(array):
        return tuple([ 'D_%010d'%d for d in array.shape ])

    #create missing dimensions
    for v in variables_save:
        for dname,d in zip( dimensions(v),v.shape ):
            if dname not in fout.dimensions.keys():
                fout.createDimension(dname, d)

    #copy variables
    for v_name, varin in fin.variables.iteritems():
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
    #parser.add_argument('files', nargs='+', type=argparse.FileType('r'))
    parser.add_argument('-a','--add', nargs='+', type=argparse.FileType('r'), help="Add the real part to the final db")
    parser.add_argument('-s','--substract', nargs='+', type=argparse.FileType('r'), help="Substract the real part to the final db", default=[])
    parser.add_argument('-ai','--addimg', nargs='+', type=argparse.FileType('r'), help="Add the imaginary part to the final db",default=[])
    parser.add_argument('-o','--output',                       help='Output filename', default='ndb_out.QP')
    parser.add_argument('-v','--verbose', action="store_true", help='Verbose mode')
    args = parser.parse_args()

#    if args.files is None:
#        parser.print_help()
#        exit()

    output  = args.output
    add     = args.add
    substract = args.substract
    addimg  = args.addimg
    #files   = args.files
    verbose = args.verbose
    add_qp(output,add,substract,addimg,verbose)




