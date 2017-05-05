# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yambopy.netcdf import *
from cmath import polar 
from yambopy.units import *
from itertools import product

class YamboExcitonDB(YamboSaveDB):
    """ Read the excitonic states database from yambo
    """
    def __init__(self,lattice,filename='ndb.BS_diago_Q01',path='.'):
        self.lattice = lattice
        self.filename = filename
        self.path = path
        self.get_database()
        
    def get_database(self):
        """ Load the diago database to memory
        """
        try:
            filename = "%s/%s"%(self.path,self.filename)
            db = Dataset(filename)
        except:
            print "failed to read database %s"%filename
            exit(1)
        if 'BS_left_Residuals' in db.variables.keys():
            #residuals
            rel,iml = db['BS_left_Residuals'][:].T
            rer,imr = db['BS_right_Residuals'][:].T
            self.l_residual = rel+iml*I
            self.r_residual = rer+imr*I
        if 'BS_Residuals' in db.variables.keys():
            #residuals
            rel,iml,rer,imr = db['BS_Residuals'][:].T
            self.l_residual = rel+iml*I
            self.r_residual = rer+imr*I
        #energies
        eig =  db['BS_Energies'][:]*ha2ev
        self.eigenvalues = eig[:,0]+eig[:,1]*I
        #eigenvectors
        eiv = db['BS_EIGENSTATES'][:]
        eiv = eiv[:,:,0] + eiv[:,:,1]*I
        self.eigenvectors = eiv
        #indexes
        self.table = db['BS_TABLE'][:].T.astype(int)

        #transitions dictionary
        #bs table k, v, c
        self.unique_vbands = np.unique(self.table[:,1]-1)
        self.unique_cbands = np.unique(self.table[:,2]-1)

        #initialize empty dictionary
        transitions_v_to_c = dict([ ((v,c),[]) for v,c in product(self.unique_vbands,self.unique_cbands) ])

        #add elements to dictionary
        for eh,kvc in enumerate(self.table-1):
            k,v,c = kvc 
            transitions_v_to_c[(v,c)].append((k,eh))

        #make an array 
        for t,v in transitions_v_to_c.items():
            transitions_v_to_c[t] = np.array(v)

        self.transitions_v_to_c = transitions_v_to_c 
        self.nexcitons    = len(self.eigenvalues)
        self.ntransitions = len(self.table)
        db.close()
    
    def exciton_bs(self,energies,path,excitons=(0,),debug=False):
        """
        Calculate exciton band-structure
            
            Arguments:
            energies -> can be an instance of YamboSaveDB or YamboQBDB
            path     -> path in reduced coordinates in which to plot the band structure
            exciton  -> exciton index to plot
        """
        if isinstance(excitons, int):
            excitons = (excitons,)
        #get full kmesh
        kpoints = self.lattice.red_kpoints
        path = np.array(path)

        kpoints_rep, kpoints_idx_rep = replicate_red_kmesh(kpoints,repx=range(-1,2),repy=range(-1,2),repz=range(-1,2))
        band_indexes = get_path(kpoints_rep,path)
        band_kpoints = kpoints_rep[band_indexes] 
        band_indexes = kpoints_idx_rep[band_indexes]

        if debug:
            for i,k in zip(band_indexes,band_kpoints):
                x,y,z = k
                plt.text(x,y,i) 
            plt.scatter(kpoints_rep[:,0],kpoints_rep[:,1])
            plt.plot(path[:,0],path[:,1],c='r')
            plt.scatter(band_kpoints[:,0],band_kpoints[:,1])
            plt.show()
            exit()

        #get eigenvalues along the path
        if isinstance(energies,YamboSaveDB):
            #exapnd eigenvalues to the bull brillouin zone
            energies = energies.eigenvalues[self.lattice.kpoints_indexes]

        #get weight of state in each band
        weights = np.zeros(energies.shape)
        for exciton in excitons:
            #get the eigenstate
            eivec = self.eigenvectors[exciton]

            for t,transitions in self.transitions_v_to_c.items():
                c,v = t
                iks, ehs = transitions.T
                weights[iks,c] += np.abs(eivec[ehs])
                weights[iks,v] += np.abs(eivec[ehs])

        energies = energies[band_indexes]
        weights  = weights[band_indexes]
        
        return np.array(band_kpoints), energies, weights 

    def plot_exciton_bs(self,ax,energies,path,excitons,size=500,space='bands'):
        """
        Plot the excitons
        
            Arguments:
            ax -> axis extance of matplotlib to add the plot to
        """
        bands_kpoints, energies, weights = self.exciton_bs(energies, path, excitons)
        
        weights /= np.max(weights)

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk-1]-bands_kpoints[nk])
            bands_distances.append(distance)

        for v,c in product(self.unique_vbands,self.unique_cbands):
            if space=='bands':
                ax.plot(bands_distances, energies[:,c], c='b')
                ax.plot(bands_distances, energies[:,v], c='b')
                ax.scatter(bands_distances, energies[:,c], s=weights[:,c]*size, c='r')
                ax.scatter(bands_distances, energies[:,v], s=weights[:,v]*size, c='r')
            else:
                ax.plot(bands_distances, energies[:,c]-energies[:,v], c='b')
                ax.scatter(bands_distances, energies[:,c]-energies[:,v], s=weights[:,c]*size, c='r')

    def get_amplitudes_phases(self,excitons=(0,)):
        """ get the excitonic amplitudes and phases
        """
        if isinstance(excitons, int):
            excitons = (excitons,)
        car_kpoints = self.lattice.car_kpoints
       
        nkpoints = len(car_kpoints)
        amplitudes = np.zeros([nkpoints])
        phases     = np.zeros([nkpoints],dtype=np.complex64)
        for exciton in excitons:
            #the the eigenstate
            eivec = self.eigenvectors[exciton]
           
            total = 0
            for eh,kvc in enumerate(self.table):
                ikbz, v, c = kvc-1
                Acvk = eivec[eh]
                phases[ikbz]     += Acvk
                amplitudes[ikbz] += np.abs(Acvk)

        return car_kpoints, amplitudes, np.angle(phases)

    def chi(self,dipoles,dir=0,emin=0,emax=8,estep=0.01,broad=0.1,nexcitons='all'):
        """
        Calculate the dielectric response function using excitonic states
        """
        if nexcitons == 'all': nexcitons = self.nexcitons

        #energy range
        w = np.arange(emin,emax,estep,dtype=np.float32)
        nenergies = len(w)
        
        print "energy range: %lf -> +%lf -> %lf "%(emin,estep,emax)
        print "energy steps: %lf"%nenergies

        #initialize the susceptibility intensity
        chi = np.zeros([len(w)],dtype=np.complex64)

        #calculate exciton-light coupling
        print "calculate exciton-light coupling"
        EL1,EL2 = self.project1(dipoles.dipoles[:,dir],nexcitons) 

        #get dipole
        #dip1 = self.l_residual
        #dip2 = self.r_residual

        #iterate over the excitonic states
        for s in xrange(nexcitons):
            #get exciton energy
            es = self.eigenvalues[s]
 
            #calculate the green's functions
            G1 = 1/(   w - es - broad*I) 
            G2 = 1/( - w - es - broad*I)

            r = EL1[s]*EL2[s]
            chi += r*G1 + r*G2

        return w,chi

    def __str__(self):
        s  = "number of excitons: %d\n"%self.nexcitons
        s += "number of excitons: %d\n"%self.ntransitions
        return s
