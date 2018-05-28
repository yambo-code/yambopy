from yambopy import *
import numpy as np
import shutil
import os
from netCDF4 import Dataset

def abs2(x):
    return x.real**2 + x.imag**2

class YamboWFDB():
    def __init__(self,savedb,path=None,save='SAVE',filename='ns.wf'):
        """
        load wavefunction from yambo
        """
        if path is None:
            self.path = save
        else:
            self.path = path+'/SAVE'
        self.filename = filename
        
        #take some data from savedb
        self.savedb   = savedb
        self.wfcgrid  = savedb.wfcgrid
        self.gvectors = savedb.gvectors
        self.kpoints  = savedb.kpts_car
        self.lat  = savedb.lat
        self.rlat = savedb.rlat
       
        #read wf 
        self.read()
        self.nkpoints, self.nspin, self.ng, self.nbands = self.wf.shape

    def read(self):
        path = self.path
        filename = self.filename

        wf = []
        nk = 1
        while True:
            try:
                fname = "%s/%s_fragments_%d_1"%(path,filename,nk)
                database = Dataset(fname)
                re = database.variables['WF_REAL_COMPONENTS_@_K%d_BAND_GRP_1'%nk][:]
                im = database.variables['WF_IM_COMPONENTS_@_K%d_BAND_GRP_1'%nk][:]
                a = re+1j*im
                wf.append(a)
                nk+=1
            except:
                if nk==1:
                    raise IOError('Could not read %s'%fname)
                break
        self.wf = np.array(wf)
        self.nkpoints, self.nspin, self.ng, self.nbands = self.wf.shape

    def get_wf_gvecs(self,kpoint=0):
        """
        Get the indexes of teh wavefunctions
        """

        #create array for fft
        indexes = self.wfcgrid[kpoint]
        indexes = indexes[indexes > 0] #remove componnents that do not belong 
        gvecs = self.gvectors[indexes]

        return gvecs

    def write(self,path):
        """
        write the wavefunctions in new files
        """
        if os.path.isdir(path): shutil.rmtree(path)
        os.mkdir(path)

        #copy all the files
        oldpath = self.path
        filename = self.filename
        shutil.copyfile("%s/%s"%(oldpath,filename),"%s/%s"%(path,filename))
        for nk in xrange(self.nkpoints):
            fname = "%s_fragments_%d_1"%(filename,nk+1)
            shutil.copyfile("%s/%s"%(oldpath,fname),"%s/%s"%(path,fname))

        #edit with the new wfs
        wf = self.wf
        for nk in xrange(self.nkpoints):
            fname = "%s_fragments_%d_1"%(filename,nk+1)
            database = Dataset("%s/%s"%(path,fname),'r+')
            database.variables['WF_REAL_COMPONENTS_@_K%d_BAND_GRP_1'%(nk+1)][:] = wf[nk].real
            database.variables['WF_IM_COMPONENTS_@_K%d_BAND_GRP_1'%(nk+1)][:] = wf[nk].imag
            database.close()
        print 'new wavefunctions written in %s'%path

    def __str__(self):
        s = ""
        s += "nkpoints: %4d\n"%self.nkpoints
        s += "nspin:    %4d\n"%self.nspin
        s += "nbands:   %4d\n"%self.nbands
        s += "ng:       %4d\n"%self.ng
        return s


if __name__ == "__main__":
    ywf = YamboWFDB(path='database')
