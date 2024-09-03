#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC
#
# This file is part of the yambopy project
#
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from yambopy.lattice import red_car
from yambopy.kpoints import get_path
from yambopy.dbs.latticedb import YamboLatticeDB
from itertools import product
from netCDF4 import Dataset

class YamboExcitonWeight(object):
    """
    Class to read the excitonic weight writen by ypp
    """

    def __init__(self,filename,save='SAVE',path='.'):
        #read save database

        # ! I have commented line 20 because is was broken. Pls check

        ysave = YamboLatticeDB.from_db_file(folder=save,filename='ns.db1')
        self.sym_car  = ysave.sym_car
        self.kpts_car = ysave.car_kpoints #kpts_car
        self.rlat     = ysave.rlat
        self.lat = ysave.lat

        #super(YamboLatticeDB,self).from_db_file(folder=save,filename='ns.db1')

        #read excitons file
        self.excitons = np.loadtxt(filename)

        self.weights = None

    def write_irr(self,filename="irr.dat"):
        """ write the list of kpoints from the irreducible brillouin zone
        """
        f = open("irr.dat",'w')
        for k in self.kpts_car:
            f.write(("%12.8lf"*3)%(k[1],k[0],k[2])+"\n")
        f.close()

    def write_full(self,filename="full.dat"):
        """ write the list of kpoints in the full brillouin zone
        """
        #generate all the possible points
        kpoints = self.kpts_car
        kpts = []
        for k in self.kpts_car:
            for sym in self.sym_car:
                kpts.append(np.dot(sym,k))

        f = open("full.dat",'w')
        for q in kpts:
            f.write(("%12.8lf "*3)%tuple(q)+"\n")
        f.close()

    def get_data(self):
        qpts, weights, transitions = self.calc_kpts_weights()
        return { "qpts": qpts,
                 "weights": weights,
                 "lattice": self.lat,
                 "reciprocal_lattice": self.rlat,
                 "transitions": transitions
                 }

    def calc_kpts_weights(self,repx=list(range(-1,2)),repy=list(range(-1,2)),repz=list(range(-1,2))):
        """ Calculate the weights and kpoints of the excitons
        """
        self.weights     = dict()

        #first run set everything to zero
        for line in self.excitons:
            v,c,k,sym,w,e = line
            self.weights[(int(k),int(sym))] = 0

        #add weights
        for line in self.excitons:
            v,c,k,sym,w,e = line
            self.weights[(int(k),int(sym))] += w

        #rename symmetries and kpoints
        sym = self.sym_car
        kpoints = self.kpts_car

        qpts     = []
        kidx     = []
        weights  = []

        for r in product(repx,repy,repz):
          for k,s in list(self.weights.keys()):
            w   = self.weights[(k,s)]
            weights.append( w )
            qpt = np.dot(sym[s-1],kpoints[k-1])+red_car([r],self.rlat)[0]
            qpts.append( qpt )
            kidx.append( k )

        return np.array(qpts), np.array(weights)

    def calc_kpts_transitions(self,repx=list(range(-1,2)),repy=list(range(-1,2)),repz=list(range(-1,2)),debug=False):
        """ Calculate the transitions and kpoints of the excitons
        """
        self.weights     = dict()
        self.transitions = dict()
        self.transitions_v_to_c = dict()

        #first run set everything to zero
        for line in self.excitons:
            v,c,k,sym,w,e = line
            self.weights[(int(k),int(sym))] = 0
            self.transitions[(int(v),int(c),int(k),int(sym))] = 0
            self.transitions_v_to_c[(int(v),int(c))] = 0

        #add weights
        for line in self.excitons:
            v,c,k,sym,w,e = line
            self.weights[(int(k),int(sym))] += w

        #add transitions
        for line in self.excitons:
            v,c,k,sym,w,e = line
            self.transitions[(int(v),int(c),int(k),int(sym))] += w

        #add percentage of a given v => c transition
        norm = sum(self.excitons[:,4])
        for v,c,k,s in list(self.transitions.keys()):
          self.transitions_v_to_c[(int(v),int(c))] += self.transitions[(v,c,k,s)]
        if debug: print('transitions (valence > condution):')
        for v,c in self.transitions_v_to_c:
          self.transitions_v_to_c[(v,c)] = self.transitions_v_to_c[(v,c)]/norm
          if debug: print('%3d > %3d'%(v,c))

        #rename symmetries and kpoints
        sym = self.sym_car
        kpoints = self.kpts_car

        qpts     = []
        kidx     = []
        t_v_c    = []

        for r in product(repx,repy,repz):
          for k,s in list(self.weights.keys()):
            qpt = np.dot(sym[s-1],kpoints[k-1])+red_car([r],self.rlat)[0]
            qpts.append( qpt )
            kidx.append( k )
            #print (v_ref,c_ref,k,s)
            #aux.append(self.transitions[(v_ref,c_ref,k,s)])

        for v_ref,c_ref in list(self.transitions_v_to_c.keys()):
          aux = []
          for r in product(repx,repy,repz):
            for k,s in list(self.weights.keys()):
              aux.append(self.transitions[(v_ref,c_ref,k,s)])
          t_v_c.append(np.array(aux))

        return np.array(qpts), t_v_c, np.array(kidx)

    def plot_contour(self,resX=500,resY=500):
        """ plot a contour
            idea taken from http://stackoverflow.com/questions/18764814/make-contour-of-scatter
        """
        kpts, z = self.calc_kpts_weights()
        x,y = kpts[:,0],kpts[:,1]
        xi = np.linspace(min(x), max(x), resX)
        yi = np.linspace(min(y), max(y), resY)
        Z = griddata(x, y, z, xi, yi, interp='cubic')
        X, Y = np.meshgrid(xi, yi)

        plt.contourf(X, Y, Z, cmap='gist_heat_r')
        plt.show()

    def plot_weights(self,ax,size=20,marker='H',set_origin=0.0,lim=0.2,cmap='viridis',log_scale=False,set_maximum=1.0):
        """

        Plot the weights in a scatter plot of this exciton (1st version tuned by A. Molina-Sanchez)
        Options:
        cmap : colormap. Default viridis 
        log_scale : Logarithmic scale for the intensity (True or False)
        set_maximum : Only applied for linear scale. Apply a cut for a selected intensity (values between 0 and 1)
        Further development: Option for the colorbar
 
        """
        from numpy import sqrt
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors

        """
        These options can be decided by the user. In this first version we just:
        remove axis
        """
        ax.set_aspect('equal')
        ax.set_xlim(-lim,lim)
        ax.set_ylim(-lim,lim)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

        
        kpts, weights = self.calc_kpts_weights()

        if log_scale == True:
           norm = colors.LogNorm(vmin=weights.min(),vmax=weights.max())
        else:
           if abs(set_maximum)>1.:
              set_maximum = 1.
           norm = colors.Normalize(vmin=weights.min(),vmax=abs(set_maximum)*weights.max())

        cmap = plt.get_cmap(cmap)

        ax.scatter(kpts[:,0]-set_origin, kpts[:,1]-set_origin, s=size, marker=marker, color=cmap(norm(weights)))

    def __str__(self):
        s = ""
        s += "reciprocal lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.rlat])+"\n"
        s += "lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.lat])+"\n"
        s += "alat:\n"
        s += ("%12.8lf "*3)%tuple(self.alat)+"\n"
        return s

    def plot_transitions(self,size=30,lim=0.2):
        """
        Plot the weight of a given transition in a scatter plot of this exciton.
        My idea is to associate for each transition a color and to plot in a different plot
        """

        from numpy import sqrt
        cmap = plt.get_cmap("gist_heat_r")

        fig = plt.figure(figsize=(10,10))
        kpts, t_v_c, _ = self.calc_kpts_transitions()
        for individual in t_v_c:
          plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', color=[cmap(sqrt(c)) for c in individual])

        plt.xlim([-lim,lim])
        plt.ylim([-lim,lim])
        ax = plt.axes()
        ax.set_aspect('equal')
        plt.show()

    def plot_exciton_bs(self,ax,path,nbands='all',space='transition',color='#1f77b4'):
        """
        Plot the excitonic weights of a given transition in the band-structure
        """
        kpts, t_v_c, kidx = self.calc_kpts_transitions(repx=list(range(1)),repy=list(range(1)),repz=list(range(1)))
        t_v_c = np.array(t_v_c)

        #get_path is provided by savedb
        bands_kpoints, bands_indexes, path_car = get_path(kpts,self.rlat,None,path)

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk-1]-bands_kpoints[nk])
            bands_distances.append(distance)

        #get energies at these k-points
        eig = self.eigenvalues[kidx-1]
        eig = eig[bands_indexes]
        transition_weight  = t_v_c[:,bands_indexes]
        if nbands == 'all': nbands = self.nbands
        for tw,t in zip(transition_weight,list(self.transitions_v_to_c.keys())):
            v,c = t
            if space == 'transition':
                ax.plot(bands_distances,eig[:,c-1]-eig[:,v-1])
                ax.scatter(bands_distances,eig[:,c-1]-eig[:,v-1],s=tw*1e4)
            else:
                ax.plot(bands_distances,eig[:,c-1],c=color)
                ax.plot(bands_distances,eig[:,v-1],c=color)
                ax.scatter(bands_distances,eig[:,c-1],s=tw*1e4,c=color)
                ax.scatter(bands_distances,eig[:,v-1],s=tw*1e4,c=color)

    def __str__(self):
        s = ""
        s += "reciprocal lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.rlat])+"\n"
        s += "lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.lat])+"\n"
        s += "alat:\n"
        s += ("%12.8lf "*3)%tuple(self.alat)+"\n"
        return s
