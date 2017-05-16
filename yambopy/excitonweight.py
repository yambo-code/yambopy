from __future__ import print_function
# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from builtins import range
from yambopy import *
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
        qpts, weights = self.calc_kpts_weights()
        return { "qpts": qpts,
                 "weights": weights,
                 "lattice": self.lat,
                 "reciprocal_lattice": self.rlat }

    def calc_kpts_weights(self,repx=list(range(3)),repy=list(range(3)),repz=list(range(3))):
        """ Calculate the weights and kpoints of the excitons
        """
        self.weights = dict()

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

        qpts = []
        weights = []
        for r in product(repx,repy,repz):
            for k,s in list(self.weights.keys()):
                w   = self.weights[(k,s)]
                weights.append( w )
                qpt = np.dot(sym[s-1],kpoints[k-1])+red_car([r],self.rlat)[0]
                qpts.append( qpt )
        return np.array(qpts), np.array(weights)

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

    def plot_weights(self,size=20):
        """ Plot the weights in a scatter plot of this exciton
        """
        cmap = plt.get_cmap("gist_heat_r")

        fig = plt.figure(figsize=(20,20))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
        kpts, weights = self.calc_kpts_weights()
        plt.scatter(kpts[:,0], kpts[:,1], s=size, marker='H', color=[cmap(sqrt(c)) for c in weights])
        plt.axes().set_aspect('equal', 'datalim')
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

if __name__ == "__main__":
    ye = YamboExciton('o-yambo.exc_weights_at_1_02')
    print(ye)
    ye.write_irr()
    ye.write_full()
    #ye.plot_contour()
    ye.plot_weights()
