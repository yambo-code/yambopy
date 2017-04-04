# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.netcdf import *
from itertools import product

class YamboExcitonWeight(YamboSaveDB):
    def __init__(self,filename,save='SAVE',path='.'):
        #read save database
        YamboSaveDB.__init__(self,save=save)

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

    def calc_kpts_weights(self,repx=range(-1,2),repy=range(-1,2),repz=range(-1,2)):
        """ Calculate the weights and kpoints of the excitons
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
            self.transitions[(int(v),int(c),int(k),int(sym))] += w

        #add percentage of a given v => c transition
        norm = sum(self.excitons[:,4])
        for v,c,k,s in self.transitions.keys():
          self.transitions_v_to_c[(int(v),int(c))] += self.transitions[(v,c,k,s)]
        for v,c in self.transitions_v_to_c:
          self.transitions_v_to_c[(v,c)] = self.transitions_v_to_c[(v,c)]/norm
          print('v ', v,' ==>>>', ' c ',c)

        #rename symmetries and kpoints
        sym = self.sym_car
        kpoints = self.kpts_car

        qpts     = []
        kidx     = []
        weights  = []
        t_v_c    = []

        for r in product(repx,repy,repz):
          for k,s in self.weights.keys():
            w   = self.weights[(k,s)]
            weights.append( w )
            qpt = np.dot(sym[s-1],kpoints[k-1])+red_car([r],self.rlat)[0]
            qpts.append( qpt )
            kidx.append( k )
            #print (v_ref,c_ref,k,s)
            #aux.append(self.transitions[(v_ref,c_ref,k,s)])
         
        for v_ref,c_ref in self.transitions_v_to_c.keys():
          aux = []
          for r in product(repx,repy,repz):
            for k,s in self.weights.keys():
              aux.append(self.transitions[(v_ref,c_ref,k,s)])
          t_v_c.append(np.array(aux)) 

        return np.array(qpts), np.array(weights), t_v_c, np.array(kidx)

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

    def plot_weights(self,size=30,lim=0.2):
        """
        Plot the weights in a scatter plot of this exciton
        """
        from numpy import sqrt
        cmap = plt.get_cmap("gist_heat_r")

        fig = plt.figure(figsize=(20,20))
        kpts, weights, transitions, _ = self.calc_kpts_weights()
        plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', color=[cmap(sqrt(c)) for c in weights])

        plt.xlim([-lim,lim])
        plt.ylim([-lim,lim])

        ax = plt.axes()
        ax.set_aspect('equal')
        plt.show()

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
        kpts, weights, t_v_c, _ = self.calc_kpts_weights()
        for individual in t_v_c:   
          plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', color=[cmap(sqrt(c)) for c in individual])

        plt.xlim([-lim,lim])
        plt.ylim([-lim,lim])
        ax = plt.axes()
        ax.set_aspect('equal')
        plt.show()

    def plot_exciton_bs(self,path,nbands='all'):
        """
        Plot the excitonic weights in the band-structure
        """
        kpts, weights, t_v_c, kidx = self.calc_kpts_weights(repx=xrange(1),repy=xrange(1),repz=xrange(1))
        t_v_c = np.array(t_v_c)

        #get_path is provided by savedb
        bands_kpoints, bands_indexes, path_car = self.get_path(path,kpts=kpts)

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
        for tw,t in zip(transition_weight,self.transitions_v_to_c.keys()):
            v,c = t
            plt.plot(bands_distances,eig[:,c-1]-eig[:,v-1])
            plt.scatter(bands_distances,eig[:,c-1]-eig[:,v-1],s=tw*1e4)
        plt.show()

    def __str__(self):
        s = ""
        s += "reciprocal lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.rlat])+"\n"
        s += "lattice:\n"
        s += "\n".join([("%12.8lf "*3)%tuple(r) for r in self.lat])+"\n"
        s += "alat:\n"
        s += ("%12.8lf "*3)%tuple(self.alat)+"\n"
        return s

    def plot_exciton_band(self,path,prefix_gw='gw',json_filename='gw'):
        """ For each transition gives a color and plot all of them in
        the LDA or GW band structure. For instance (v=4, c=5, red),
        (v=4,c=6, blue). Same weight associated to valence and cond.
        """
        band_gw = YamboAnalyser('gw')#prefix_gw)
        
        bands_kpoints, bands_indexes, bands_highsym_qpts = self.get_path(path,json_filename)
        print( bands_highsym_qpts ) 


if __name__ == "__main__":
    ye = YamboExciton('o-yambo.exc_weights_at_1_02')
    print ye
    ye.write_irr()
    ye.write_full()
    #ye.plot_contour()
    ye.plot_weights()
