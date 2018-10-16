# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import os
from itertools import product
from yambopy import *
from cmath import polar 
from yambopy.units import *
from yambopy.plot.plotting import add_fig_kwargs
from yambopy.lattice import replicate_red_kmesh, calculate_distances, get_path

class ExcitonList():
    """
    Container class to perform operations on lists of excitons
    """
    def __init__(self,excitonlist):
        self.excitonlist = excitonlist

    def __str__(self):
        lines = []; app = lines.append
        for exciton in self.excitonlist:
            app( str(exciton.get_string(singleline=True)) )
        return "\n".join(lines)

class Exciton():
    """
    Basic container of data for a single exciton
    TODO: classify the excitons according to their symmetry
    """
    def __init__(self,energy,intensity,degeneracy,coeffs=None,wf=None):
        self.energy = energy
        self.intensity = intensity
        self.degeneracy = degeneracy
        self.coeffs = coeffs
        self.wf = wf

    def get_string(self,singleline=False):
        lines = []; app = lines.append
        app( 'energy:     %8.4'%self.energy )
        app( 'intensity:  %8.4'%self.intensity )
        app( 'degeneracy: %8d'%self.degeneracy )
        if singleline: return "".join(lines)
        return "\n".join(lines)
    
    def __str__(self):
        return self.get_string()

class YamboExcitonDB(YamboSaveDB):
    """ Read the excitonic states database from yambo
    """
    def __init__(self,lattice,eigenvalues,l_residual,r_residual,table=None,eigenvectors=None):
        self.lattice = lattice
        self.eigenvalues = eigenvalues
        self.l_residual = l_residual
        self.r_residual = r_residual
        #optional
        self.table = table
        self.eigenvectors = eigenvectors

    @classmethod
    def from_db_file(cls,lattice,filename='ndb.BS_diago_Q01',folder='.'):
        """ initialize this class from a file
        """
        path_filename = os.path.join(folder,filename)
        if not os.path.isfile(path_filename):
            raise FileNotFoundError("File %s not found in YamboExcitonDB"%path_filename)

        with Dataset(path_filename) as database:
            if 'BS_left_Residuals' in list(database.variables.keys()):
                #residuals
                rel,iml = database.variables['BS_left_Residuals'][:].T
                rer,imr = database.variables['BS_right_Residuals'][:].T
                l_residual = rel+iml*I
                r_residual = rer+imr*I
            if 'BS_Residuals' in list(database.variables.keys()):
                #residuals
                rel,iml,rer,imr = database.variables['BS_Residuals'][:].T
                l_residual = rel+iml*I
                r_residual = rer+imr*I

            #energies
            eig =  database.variables['BS_Energies'][:]*ha2ev
            eigenvalues = eig[:,0]+eig[:,1]*I
                
            #eigenvectors
            table = None
            eigenvectors = None
            if 'BS_EIGENSTATES' in database.variables:
                eiv = database.variables['BS_EIGENSTATES'][:]
                eiv = eiv[:,:,0] + eiv[:,:,1]*I
                eigenvectors = eiv
                table = database.variables['BS_TABLE'][:].T.astype(int)

            table = table
            eigenvectors = eigenvectors

        return cls(lattice,eigenvalues,l_residual,r_residual,table=table,eigenvectors=eigenvectors)

    @property
    def unique_vbands(self):
        return np.unique(self.table[:,1]-1)

    @property
    def unique_cbands(self):
        return np.unique(self.table[:,2]-1)

    @property
    def transitions_v_to_c(self):
        """Compute transitions from valence to conduction"""
        if hasattr(self,"_transitions_v_to_c"): return self._transitions_v_to_c
        uniq_v = self.unique_vbands
        uniq_c = self.unique_cbands
        transitions_v_to_c = dict([ ((v,c),[]) for v,c in product(uniq_v,uniq_c) ])

        #add elements to dictionary
        for eh,kvc in enumerate(self.table-1):
            k,v,c = kvc 
            transitions_v_to_c[(v,c)].append((k,eh))

        #make an array 
        for t,v in list(transitions_v_to_c.items()):
            if len(np.array(v)):
                transitions_v_to_c[t] = np.array(v)
            else:
                del transitions_v_to_c[t]

        self._transitions_v_to_c = transitions_v_to_c 
        return transitions_v_to_c
 
    @property
    def ntransitions(self):
        return len(self.table)

    @property
    def nexcitons(self):
        return len(self.eigenvalues)

    def write_sorted(self,prefix='yambo'):
        """
        Write the sorted energies and intensities to a file
        """
        #get intensities
        eig = self.eigenvalues.real
        intensities = self.get_intensities()

        #get sorted energies
        sort_e, sort_i = self.get_sorted()     

        #write excitons sorted by energy
        with open('%s_E.dat'%prefix, 'w') as f:
            for e,n in sort_e:
                f.write("%3d %12.8lf %12.8e\n"%(n+1,e,intensities[n])) 

        #write excitons sorted by intensities
        with open('%s_I.dat'%prefix,'w') as f:
            for i,n in sort_i:
                f.write("%3d %12.8lf %12.8e\n"%(n+1,eig[n],i)) 

    def get_nondegenerate(self,eps=1e-4):
        """
        get a list of non-degenerate excitons
        """
        non_deg_e   = [0]
        non_deg_idx = [] 

        #iterate over the energies
        for n,e in enumerate(self.eigenvalues):
            if not np.isclose(e,non_deg_e[-1],atol=eps):
                non_deg_e.append(e)
                non_deg_idx.append(n)

        return np.array(non_deg_e[1:]), np.array(non_deg_idx)

    def get_intensities(self):
        """
        get the intensities of the excitons
        """
        intensities = abs2(self.l_residual*self.r_residual)
        intensities /= np.max(intensities)
        return intensities

    def get_sorted(self):
        """
        Return the excitonic weights sorted according to energy and intensity
        """
        #get intensities
        eig = self.eigenvalues.real
        intensities = self.get_intensities()

        #list ordered with energy
        sort_e = sorted(zip(eig, list(range(self.nexcitons))))

        #list ordered with intensity
        sort_i = sorted(zip(intensities, list(range(self.nexcitons))),reverse=True)

        return sort_e, sort_i 

    def get_degenerate(self,index,eps=1e-4):
        """
        Get degenerate excitons
        """
        energy = self.eigenvalues[index-1]
        excitons = [] 
        for n,e in enumerate(self.eigenvalues):
            if np.isclose(energy,e,atol=eps):
                excitons.append(n+1)
        return excitons

    def exciton_bs(self,energies,path,excitons=(0,),debug=False):
        """
        Calculate exciton band-structure
            
            Arguments:
            energies -> can be an instance of YamboSaveDB or YamboQBDB
            path     -> path in reduced coordinates in which to plot the band structure
            exciton  -> exciton index to plot
        """
        if self.eigenvectors is None:
            raise ValueError('This database does not contain Excitonic states,'
                              'please re-run the yambo BSE calculation with the WRbsWF option in the input file.')
        if isinstance(excitons, int):
            excitons = (excitons,)
        #get full kmesh
        kpoints = self.lattice.red_kpoints
        path = np.array(path)

        kpoints_rep, kpoints_idx_rep = replicate_red_kmesh(kpoints,repx=list(range(-1,2)),repy=list(range(-1,2)),repz=list(range(-1,2)))
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
        if isinstance(energies,(YamboSaveDB,YamboElectronsDB)):
            #expand eigenvalues to the bull brillouin zone
            energies = energies.eigenvalues[self.lattice.kpoints_indexes]

        elif isinstance(energies,YamboQPDB):
            #expand the quasiparticle energies to the bull brillouin zone
            energies = energies.eigenvalues_qp[self.lattice.kpoints_indexes]
        else:
            raise ValueError("Energies argument must be an instance of YamboSaveDB,"
                             "YamboElectronsDB or YamboQPDB. Got %s"%(type(energies)))

        #get weight of state in each band
        weights = np.zeros(energies.shape)
        for exciton in excitons:
            #get the eigenstate
            eivec = self.eigenvectors[exciton-1]

            #add weights
            for t,transitions in list(self.transitions_v_to_c.items()):
                c,v = t
                iks, ehs = transitions.T
                weights[iks,c] += abs2(eivec[ehs])
                weights[iks,v] += abs2(eivec[ehs])

        energies = energies[band_indexes]
        weights  = weights[band_indexes]

        #make top valence band to be zero
        energies -= max(energies[:,max(self.unique_vbands)])
        
        return np.array(band_kpoints), energies, weights 

    def plot_exciton_bs_ax(self,ax,energies_db,path,excitons,size=500,space='bands',
                           args_scatter={'c':'b'},args_plot={'c':'r'},debug=False):
        """
        Plot the exciton band-structure
        
            Arguments:
            ax          -> axis extance of matplotlib to add the plot to
            lattice     -> Lattice database
            energies_db -> Energies database, can be either a SaveDB or QPDB
            path        -> Path in the brillouin zone
        """
        bands_kpoints, energies, weights = self.exciton_bs(energies_db, path, excitons, debug)
        
        weights /= np.max(weights)

        #calculate distances
        bands_distances = [0]
        distance = 0
        for nk in range(1,len(bands_kpoints)):
            distance += np.linalg.norm(bands_kpoints[nk-1]-bands_kpoints[nk])
            bands_distances.append(distance)

        for v,c in product(self.unique_vbands,self.unique_cbands):
            if space=='bands':
                ax.plot(bands_distances, energies[:,c], **args_plot)
                ax.plot(bands_distances, energies[:,v], **args_plot)
                ax.scatter(bands_distances, energies[:,c], s=weights[:,c]*size, **args_scatter)
                ax.scatter(bands_distances, energies[:,v], s=weights[:,v]*size, **args_scatter)
            else:
                ax.plot(bands_distances, energies[:,c]-energies[:,v], c='b')
                ax.scatter(bands_distances, energies[:,c]-energies[:,v], s=weights[:,c]*size, c='r')

        #add high-symmetry k-points vertical bars
        kpath_car = red_car(path,self.lattice.rlat)
        #calculate distances for high-symmetry points
        kpath_distances = calculate_distances( path )
        for d in kpath_distances:
            ax.axvline(d,c='k')

        ax.set
        ax.set_xlim([min(bands_distances),max(bands_distances)])
        ax.set_title("exciton %d-%d"%(excitons[0],excitons[-1]))
        return kpath_distances

    @add_fig_kwargs
    def plot_exciton_bs(self,energies_db,path,excitons,size=500,space='bands',
                        args_scatter={'c':'b'},args_plot={'c':'r'},debug=False):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_exciton_bs_ax(ax,energies_db,path,excitons,size=size,space=space,
                                args_scatter=args_scatter,args_plot=args_plot,debug=debug)
        return fig
        
    def get_amplitudes_phases(self,excitons=(0,),repx=list(range(1)),repy=list(range(1)),repz=list(range(1))):
        """ get the excitonic amplitudes and phases
        """
        if self.eigenvectors is None:
            raise ValueError('This database does not contain Excitonic states,'
                             'please re-run the yambo BSE calculation with the WRbsWF option in the input file.')
        if isinstance(excitons, int):
            excitons = (excitons,)
       
        car_kpoints = self.lattice.car_kpoints
        nkpoints = len(car_kpoints)
        print(nkpoints)
        amplitudes = np.zeros([nkpoints])
        phases     = np.zeros([nkpoints],dtype=np.complex64)
        for exciton in excitons:
            #the the eigenstate
            eivec = self.eigenvectors[exciton-1]
           
            total = 0
            for eh,kvc in enumerate(self.table):
                ikbz, v, c = kvc-1
                Acvk = eivec[eh]
                phases[ikbz]     += Acvk
                amplitudes[ikbz] += np.abs(Acvk)

        #replicate kmesh
        red_kmesh,kindx = replicate_red_kmesh(self.lattice.red_kpoints,repx=repx,repy=repy,repz=repz)
        car_kpoints = red_car(red_kmesh,self.lattice.rlat)

        return car_kpoints, amplitudes[kindx], np.angle(phases)[kindx]

    def get_chi(self,dipoles=None,dir=0,emin=0,emax=10,estep=0.01,broad=0.1,q0norm=1e-5,
            nexcitons='all',spin_degen=2,verbose=0):
        """
        Calculate the dielectric response function using excitonic states
        """
        if nexcitons == 'all': nexcitons = self.nexcitons

        #energy range
        w = np.arange(emin,emax,estep,dtype=np.float32)
        nenergies = len(w)
        
        if verbose:
            print("energy range: %lf -> +%lf -> %lf "%(emin,estep,emax))
            print("energy steps: %lf"%nenergies)

        #initialize the susceptibility intensity
        chi = np.zeros([len(w)],dtype=np.complex64)

        if dipoles is None:
            #get dipole
            EL1 = self.l_residual
            EL2 = self.r_residual
        else:
            #calculate exciton-light coupling
            if verbose: print("calculate exciton-light coupling")
            EL1,EL2 = self.project1(dipoles.dipoles[:,dir],nexcitons) 

        #iterate over the excitonic states
        for s in range(nexcitons):
            #get exciton energy
            es = self.eigenvalues[s]
 
            #calculate the green's functions
            G1 = -1/(   w - es + broad*I)
            G2 = -1/( - w - es - broad*I)

            r = EL1[s]*EL2[s]
            chi += r*G1 + r*G2

        #multiply facto
        d3k_factor = self.lattice.rlat_vol/self.lattice.nkpoints
        cofactor = spin_degen/(2*np.pi)**3 * d3k_factor * (4*np.pi)
        chi = chi*cofactor/q0norm**2

        return w,chi

    def plot_chi_ax(self,ax,reim='im',n_brightest=5,**kwargs):
        """Plot chi on a matplotlib axes"""
        w,chi = self.get_chi(**kwargs)
        #cleanup kwargs variables
        for var in ['dipoles','dir','emin','emax','estep','broad','q0norm','nexcitons','spin_degen','verbose']: 
            kwargs.pop(var,None)
        if 're' in reim: ax.plot(w,chi.real,**kwargs)
        if 'im' in reim: ax.plot(w,chi.imag,**kwargs)
        ax.set_ylabel('$Im(\chi(\omega))$')
        ax.set_xlabel('Energy (eV)')
        #plot vertical bar on the brightest excitons
        exc_e,exc_i = self.get_sorted()
        for i,idx in exc_i[:n_brightest]:
            exciton_energy,idx = exc_e[idx]
            ax.axvline(exciton_energy,c='k')
            ax.text(exciton_energy,0.1,idx,rotation=90)

    @add_fig_kwargs
    def plot_chi(self,n_brightest=5,**kwargs):
        """Produce a figure with chi"""
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot_chi_ax(ax,n_brightest=n_brightest,**kwargs)
        return fig

    def get_string(self,mark="="):
        lines = []; app = lines.append
        app( marquee(self.__class__.__name__,mark=mark) )
        app( "number of excitons:    %d"%self.nexcitons )
        if self.table is not None: app( "number of transitions: %d"%self.ntransitions )
        return '\n'.join(lines)
    
    def __str__(self):
        return self.get_string()
