#
# Merge the dielectric functions
#
import subprocess
import os
from netCDF4 import Dataset
files = subprocess.check_output('ls */ndb.em1?_*',shell=True).splitlines()

os.system('mkdir yambo')
os.system('cp 1/ndb.em1* yambo')
for filename in files:
    q,fname = filename.split('/')
    print q,fname[:-1]+q

    #open dataset
    ncin = Dataset(filename,'r')

    #new dataset
    ncout = Dataset('yambo/%s'%(fname[:-1]+q),'w')

    #copy dimensions
    for dname, the_dim in ncin.dimensions.iteritems():
        ncout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

    #copy variables
    for v_name, varin in ncin.variables.iteritems():
        v_name = v_name[:-1]+q
        outVar = ncout.createVariable(v_name, varin.datatype, varin.dimensions)
        
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        
        outVar[:] = varin[:]
        if 'FREQ_PARS_sec_iq' in v_name:
            outVar[0] = float(q)

    ncin.close()
    ncout.close()
