# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
# Unfolding of the electronic structure.
# Program adapted for reading of Quantum Espresso.
#
# Authors: Alejandro Molina-Sanchez and Henrique Miranda
# Revision from the first version of 24 of February of 2014
#
#
# This class need major revision. So far is uneficcient and a waste of resources.
# Future developments will be:
# 
# - Reading netcdf files from Yambo (faster than converting to xml)
# - Assignment of the correspondence G-vectors to g-vectors without reading the PC g-vectors
#

import xml.etree.ElementTree as ET
from qepy.lattice import Path, calculate_distances
from qepy.pwxml import PwXML
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from numpy import array, sqrt, cross, dot, arange, zeros, around
from sys import stdout

HatoeV = 27.2107

class UnfoldingYambo():
    """ Class to unfold the electronic structure of supercell into the
    original primitive cell. Adapted to Quantum espresso XML files.
    """

    _eig_xml  = 'eigenval.xml'
    _gkv_xml  = 'gkvectors.xml'
    _gkv_dat  = 'gkvectors.dat'
    _evc_xml  = 'evc.xml'
    _evc_dat  = 'evc.dat'
    _evc1_xml = 'evc1.xml'
    _evc2_xml = 'evc2.xml'
    _evc1_dat = 'evc1.dat'
    _evc2_dat = 'evc2.dat'

    def __init__(self,prefix_pc,prefix_sc,path_pc='.',path_sc='.',verbose=0,spin="none",convert_to_xml=True,band_min=0):
        """ 
        Initialize the structure with the paths where the datafile.xml
        of the primitive and supercell
        """
        self.prefix_pc = prefix_pc
        self.prefix_sc = prefix_sc
        self.path_pc   = path_pc
        self.path_sc   = path_sc
        
        pc_xml = PwXML(prefix=self.prefix_pc,path=self.path_pc)
        sc_xml = PwXML(prefix=self.prefix_sc,path=self.path_sc)

        self.nkpoints_pc = pc_xml.nkpoints
        self.nkpoints_sc = sc_xml.nkpoints

        self.kpoints = sc_xml.kpoints

        self.nbands_pc = pc_xml.nbands 
        self.nbands_sc = sc_xml.nbands 
        self.band_min  = band_min

        if self.band_min > self.nbands_sc:
           raise Exception("Minimum of bands larger than total number of bands")

        self.cell_pc = pc_xml.cell
        self.cell_sc = sc_xml.cell

        self.rcell_pc = array(pc_xml.rcell)/pc_xml.celldm[0]   # Dont remember if we need to write the reciprocal vectors in cart. units?
        self.rcell_sc = array(sc_xml.rcell)/sc_xml.celldm[0]

        self.eigen_pc = array(pc_xml.eigen)
        self.eigen_sc = array(sc_xml.eigen)

        format_string = "%12.4lf %12.4lf %12.4lf"
        n_decs = 8
        # Angle between prim. cell and supercell

        v1 = array(pc_xml.rcell[:][0])
        v2 = array(sc_xml.rcell[:][0])

        norm_v1 = sqrt(abs(dot(v1,v1)))
        norm_v2 = sqrt(abs(dot(v2,v2)))

        cos_v1v2 = dot(v1,v2)/(norm_v1*norm_v2)
        sin_v1v2 = sqrt(dot(cross(v1,v2),cross(v1,v2)))/(norm_v1*norm_v2)

        self.rot = array([ [cos_v1v2,sin_v1v2,0.0], [-sin_v1v2,cos_v1v2,0.0], [0.0,0.0,1.0]] )

        self.projection = zeros([self.nkpoints_sc,self.nbands_sc-self.band_min])

        save_pc = YamboLatticeDB.from_db_file(folder="%s/%s.save/SAVE" % (self.path_pc,self.prefix_pc))
        save_sc = YamboLatticeDB.from_db_file(folder="%s/%s.save/SAVE" % (self.path_sc,self.prefix_sc))

        self.wf_pc = YamboWFDB(save_pc,path=self.path_pc)
        self.wf_sc = YamboWFDB(save_sc,path=self.path_sc)

        print(self.wf_pc)

        '''
    #def convert_dat_xml(self):
        if convert_to_xml == True:
           print("converting dat files to xml...")
           for ik in range(self.nkpoints_sc):
               load(ik,self.nkpoints_sc)
               #def convert_dat_xml(self):  Bring this to a function
               
               file_dat = "%s/%s.save/K%05d/%s" % (self.path_pc,self.prefix_pc,(ik + 1),self._gkv_dat)
               file_xml = "%s/%s.save/K%05d/%s" % (self.path_pc,self.prefix_pc,(ik + 1),self._gkv_xml)
               os.system('iotk convert %s %s' % (file_dat,file_xml))
   
               file_dat = "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._gkv_dat)
               file_xml = "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._gkv_xml)
               os.system('iotk convert %s %s' % (file_dat,file_xml))
               
               #file_dat = "%s/%s.save/K%05d/%s" % (self.path_pc,self.prefix_pc,(ik + 1),self._evc_dat)
               #file_xml = "%s/%s.save/K%05d/%s" % (self.path_pc,self.prefix_pc,(ik + 1),self._evc_xml)
               #os.system('iotk convert %s %s' % (file_dat,file_xml))

               if spin == "none":
            
                  file_dat = "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc_dat)
                  file_xml = "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc_xml)
                  os.system('iotk convert %s %s' % (file_dat,file_xml))
            
               if spin == "spinor":
            
                  file_dat = "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc1_dat)
                  file_xml = "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc1_xml)
                  os.system('iotk convert %s %s' % (file_dat,file_xml))
                  file_dat = "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc2_dat)
                  file_xml = "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc2_xml)
                  os.system('iotk convert %s %s' % (file_dat,file_xml))
           print("done!") 
        ''' 

         

        #gkvectors = []
        print("Dictionary of G-vectors and projection")
        for ik in range(self.nkpoints_sc):

            load(ik,self.nkpoints_sc)

            # Reading the G-vectors and g-vectors
            tree_gk_sc = ET.parse( "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._gkv_xml) )
            tree_gk_pc = ET.parse( "%s/%s.save/K%05d/%s" % (self.path_pc,self.prefix_pc,(ik + 1),self._gkv_xml) )
            root_gk_sc = tree_gk_sc.getroot()
            root_gk_pc = tree_gk_pc.getroot()

            #get the number of g-vectors
            n_gvec_sc = int(root_gk_sc.find("GRID").get('size'))
            n_gvec_pc = int(root_gk_pc.find("GRID").get('size'))
            self.ng_sc = int(n_gvec_sc/3)
            self.ng_pc = int(n_gvec_pc/3)
            g_sc= dict() 
            #print("Dimension G- and g-vectors", self.ng_sc, self.ng_pc)

            #check this point
            for GRID in root_gk_sc.findall("GRID"):
                gkold = GRID.text.split("\n")

            g_sc_int = dict()  # dictionary of integers
            #print('Reading Supercell G-vectors')
            for ig in arange(self.ng_sc):  # ATTENTION: Why was it xrange?
                #print("Reading Supercell G-vectors") 
                #load(ig,self.ng_sc)

                x,y,z = map( float, gkold[ig+1].split())
                g_sc_int[(int(x),int(y),int(z))] = ig
                w = x*self.rcell_sc[:][0] + y*self.rcell_sc[:][1] + z*self.rcell_sc[:][2] #scaling
                w = dot(self.rot,w) #rotations
                w = around(w, decimals=n_decs)+array([0,0,0]) #round and clean
                w = format_string % (w[0],w[1],w[2]) #truncation
                g_sc[w] = ig #create dictionary
    
            #print('Assigning Primitive cell g-vectors')

            g_contain = [0]*self.ng_pc

            for GRID in root_gk_pc.findall("GRID"):
                gkold = GRID.text.split("\n")

            for ig in arange(self.ng_pc):
                #load(ig,self.ng_pc)

                x,y,z = map( float, gkold[ig+1].split())
                #print(ig,int(x),int(y),int(z))
                w = x*self.rcell_pc[:][0] + y*self.rcell_pc[:][1] + z*self.rcell_pc[:][2] #scaling
                w = around(w, decimals=n_decs)+array([0,0,0]) #round and clean
                w = format_string % (w[0],w[1],w[2]) #truncation
                try:
                    g_contain[ig] = g_sc[w]
                except KeyError:
                    print("Missing k-point %d" % ig)
                    print(w)
                    #print("g_sc ")
                    #print(w,g_sc[w])
                    #print("g_contain ")
                    #print(w,g_contain[ig])
                    #g_contain[ig] = 0
            #exit()
            #print(g_contain)
            #print("ng_pc", self.ng_pc)
            #print("ng_sc", self.ng_sc)
            
            ndim_gcontain = len(g_contain)
            #print(ik)
            #print("g_contain dimension")
            #print(ndim_gcontain)
            #print("ng_pc %d" % self.ng_pc)
            #print()
            #exit()

        #evc = []
        #for ik in range(self.nkpoints_sc):

            # Reading the Super-cell Eigenvectors

            if spin == "none":

               tree_evc_sc = ET.parse( "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc_xml) )
               root_evc_sc = tree_evc_sc.getroot()
               
            if spin == "spinor":

               tree_evc1_sc = ET.parse( "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc1_xml) )
               root_evc1_sc = tree_evc1_sc.getroot()
               tree_evc2_sc = ET.parse( "%s/%s.save/K%05d/%s" % (self.path_sc,self.prefix_sc,(ik + 1),self._evc2_xml) )
               root_evc2_sc = tree_evc2_sc.getroot()

            if spin == "none":

               eivecs = []
               for ib in range(self.band_min,self.nbands_sc):
                   eivec = root_evc_sc.find("evc."+str(ib+1)).text.split("\n")
                   eivecs.append( map(lambda x: complex( float(x.split(",")[0]), float(x.split(",")[1]) ), eivec[1:-1]) )
                   if ib==0:
                      x = 0.0
                      for ig in range(self.ng_sc):
                          x += eivecs[-1][ig]*eivecs[-1][ig].conjugate()

            if spin == "spinor":

               eivecs1, eivecs2 = [], []
               for ib in range(self.band_min,self.nbands_sc):
                   #print("Reading Supercell Wave functions") 
                   #load(ib,self.nbands_sc)

                   eivec1 = root_evc1_sc.find("evc."+str(ib+1)).text.split("\n")
                   eivec2 = root_evc2_sc.find("evc."+str(ib+1)).text.split("\n")
                   eivecs1.append(list( map(lambda x: complex( float(x.split(",")[0]), float(x.split(",")[1]) ), eivec1[1:-1]) ) )
                   eivecs2.append(list( map(lambda x: complex( float(x.split(",")[0]), float(x.split(",")[1]) ), eivec2[1:-1]) ) )
               #print(eivecs1.shape) 
               #print(eivecs2.shape) 
                   #Why is this here? Was it a test?

                   #if ib==0:
                   #   x = 0.0
                   #   for ig in range(self.ng_sc):
                   #       x += eivecs1[-1][ig]*eivecs1[-1][ig].conjugate()
                   #       x += eivecs2[-1][ig]*eivecs2[-1][ig].conjugate()



                       #print(ig, abs(eivecs[-1][ig]), eivecs[-1][ig].conjugate() )
                       #if abs(eivecs[-1][ig]*eivecs[-1][ig].conjugate()) > 0.05:
                       #    print("warning",eivecs[-1][ig]*eivecs[-1][ig].conjugate())
                   #print(x)
            #evc.append(eivecs)
            #    exit()     
        #sc.convert_dat_xml()

        # Projection
            if spin == "none":

               for ib in range(self.nbands_sc-self.band_min): 
                   x = 0.0
                   for ig in range(self.ng_pc): #ndim_gcontain):
                       x += eivecs[ib][g_contain[ig]]*(eivecs[ib][g_contain[ig]].conjugate())
                    #if ib==0:
                    #   print(ib,ig,g_contain[ig],eivecs[ib][g_contain[ig]])
                    #if ib==0:
                       #print(eivecs[ib][g_contain[ig]])
                   self.projection[ik][ib] = abs(x)

            if spin == "spinor":
                  
               #print(eivecs1)
               #exit()

               for ib in range(self.nbands_sc-self.band_min): 
                   x = 0.0
                   for ig in range(self.ng_pc): #ndim_gcontain):
                       x += eivecs1[ib][g_contain[ig]]*(eivecs1[ib][g_contain[ig]].conjugate())
                       x += eivecs2[ib][g_contain[ig]]*(eivecs2[ib][g_contain[ig]].conjugate())
                    #if ib==0:
                    #   print(ib,ig,g_contain[ig],eivecs[ib][g_contain[ig]])
                    #if ib==0:
                       #print(eivecs[ib][g_contain[ig]])
                   self.projection[ik][ib] = abs(x)

        print("Done!")

        #print(self.projection)
    # Plotting adapted from PwXML (to be improved)

    def plot_eigen_ax(self,ax,path=[],xlim=(),ylim=()):

        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            ax.set_xticks( *list(zip(*path)) )
        ax.set_ylabel('E (eV)')

        #get kpoint_dists 
        kpoints_dists = calculate_distances(self.kpoints)

        #make labels
        ticks, labels = list(zip(*path))
        ax.set_xticks([kpoints_dists[t] for t in ticks])
        ax.set_xticklabels(labels)
        ax.set_ylabel('E (eV)')

        #plot vertical line
        for t in ticks:
            ax.axvline(kpoints_dists[t],c='k',lw=2)
        ax.axhline(0,c='k',lw=1)

        #plot bands
        for ib in range(self.nbands_pc):
           #ax.plot(list(range(self.nkpoints_pc)),self.eigen_pc[:,ib]*HatoeV,'k--',lw=0.5)
           ax.plot(kpoints_dists,self.eigen_pc[:,ib]*HatoeV,'k--',lw=0.5)

        for ib in range(self.nbands_sc-self.band_min):
           #ax.plot(list(range(self.nkpoints_sc)),self.eigen_sc[:,ib]*HatoeV,'r--',lw=1)
           #ax.scatter(list(range(self.nkpoints_sc)),self.eigen_sc[:,ib]*HatoeV,s=self.projection[:,ib]*20,color='r')
           ax.scatter(kpoints_dists,self.eigen_sc[:,ib+self.band_min]*HatoeV,s=self.projection[:,ib]*20,color='r')

        #plot options
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)

def load(x,n):
    bar_length = 100
    x+=1
    ratio = x/float(n)
    c = int(ratio * bar_length)
    stdout.write("["+"="*c+" "*(bar_length-c)+"] %03.3f%%" % (ratio*100))
    if (x==n): stdout.write("\n")
    stdout.flush()
    stdout.write("\r")

