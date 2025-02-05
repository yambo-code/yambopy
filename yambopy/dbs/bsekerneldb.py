# Author: Davide Romanin (revised: FP, RR)
#
# This file is part of the yambopy project
#
import os
from yambopy import *
from yambopy.units import *
from itertools import product

class YamboBSEKernelDB(object):
    """ Read the BSE Kernel database from yambo.
        It reads <t1| K |t2> where K is the kernel and t_i transition indices.
        
        Can only be used if yambo is run with parallel IO.
        
        Only supports "RESONANT" case for BSE calculation.
        TODO: support more cases
    """
    def __init__(self,lattice,kernel):
        if not isinstance(lattice,YamboLatticeDB):
            raise ValueError('Invalid type for lattice argument. It must be YamboLatticeDB')
        
        self.lattice  = lattice
        self.kernel   = kernel

    @classmethod
    def from_db_file(cls, lattice, Qpt=1, folder='.'):
        """Initialize this class from a ndb.BS_PAR_Q# file."""
        filename = 'ndb.BS_PAR_Q%d' % Qpt
        path_filename = os.path.join(folder, filename)
        if not os.path.isfile(path_filename):
            raise FileNotFoundError(f"File {path_filename} not found in YamboExcitonDB")

        with Dataset(path_filename) as database:
            if 'BSE_RESONANT' in database.variables:
                # Read as transposed since dimensions in netCDF are inverted
                reker, imker = database.variables['BSE_RESONANT'][:].T
                ker = reker + imker*I
                
                # Transform the triangular matrix to a square Hermitian matrix
                kernel = np.conjugate(np.transpose(np.triu(ker))) + np.triu(ker)
                kernel[np.diag_indices(len(kernel))] *= 0.5
                
                # Check if the kernel is Hermitian
                if not np.allclose(kernel, np.conjugate(kernel.T)):
                    raise ValueError("The constructed kernel matrix is not Hermitian")
            else:
                raise ValueError('Only BSE_RESONANT case supported so far')

        return cls(lattice, kernel)

    @property
    def ntransitions(self): return len(self.kernel)

    def consistency_BSE_BSK(self,excitons):
        """ Check that exciton and kernel dbs are consistent
        """
        if excitons.nexcitons != self.ntransitions:
            print('[WARNING] Exciton and transition spaces have different dimensions!')
        if excitons.ntransitions != self.ntransitions:
            print('[WARNING] Mismatch in ntransitions between ExcitonDB and BSEkernelDB!')        

    def get_kernel_exciton_basis(self,excitons):
        """ Switch from transition |tq>=|kc,k-qv> to excitonic |lq> basis. 
            In this basis the kernel is diagonal.
            
            <l|K|l> = sum_{t,t'}( <l|t><t|K|t'><t'|l> )
                    = sum_{t,t'}( (A^l_t)^* K_tt' A^l_t' ) 
       
            exciton: YamboExcitonDB object 
            Here t->kcv according to table from YamboExcitonDB database
        """
        kernel   = self.kernel
        Nstates  = self.ntransitions
        eivs     = excitons.eigenvectors
        self.consistency_BSE_BSK(excitons)

        # Basis transformation
        kernel_exc_basis  = np.einsum('ij,kj,ki->k', kernel, eivs, np.conj(eivs), optimize=True)
        #kernel_exc_basis = np.zeros(Nstates,dtype=complex)
        #for il in range(Nstates):
        #    kernel_exc_basis[il] = np.dot( np.conj(eivs[il]), np.dot(kernel,eivs[il]) )

        return kernel_exc_basis

    def get_kernel_value_bands(self,excitons,bands):
        """ Get value of kernel matrix elements 
            as a function of k in BZ for fixed c,v bands:
            
                K_cv(k,p) = <ck,vk-q|K|cp,vp-q>
                
            exciton: YamboExcitonDB object
            bands = [iv,ic] (NB: enumerated starting from one instead of zero) 
        """
        table  = excitons.table
        nk     = self.lattice.nkpoints
        kernel = self.kernel
        self.consistency_BSE_BSK(excitons)
        
        if bands[0] not in table[:,1] or bands[1] not in table[:,2]:
            raise ValueError('Band indices not matching available transitions')
             
        # Wcv defined on the full BZ (only a subset will be filled)
        Wcv = np.zeros((nk,nk),dtype=complex)
        # Find indices where selected valence band appears
        t_v = np.where(table[:,1]==bands[0])[0]
        # Find indices where selected conduction band appears
        t_c = np.where(table[:,2]==bands[1])[0]
        # Among those, find subset of indices where both appear together
        t_vc = [t for t in t_v if t in t_c ]

        # Iterate only on the subset
        for it1_subset, it2_subset in product(t_vc,repeat=2):
            ik = table[it1_subset][0]
            ip = table[it2_subset][0]
            Wcv[ik-1,ip-1] = kernel[it1_subset,it2_subset]
        return Wcv

    def get_kernel_value_bands_4D(self, excitons, bands_range):
        """
        Get value of kernel matrix elements as a function of (v, c, k, p)
        where v and c are valence and conduction bands, k and p are k-point indices.

        Args:
        excitons: YamboExcitonDB object containing exciton data and the lattice structure
        bands_range: Tuple (min_band, max_band) defining the range of band indices

        Returns:
        W: A 4D numpy array with dimensions corresponding to (v, c, k, p)
        """
        table = excitons.table
        nk = self.lattice.nkpoints
        kernel = self.kernel
        self.consistency_BSE_BSK(excitons)

        min_band, max_band = bands_range
        nbands = max_band

        # Initialize the 4D array for kernel values with dimensions covering the band range
        W = np.zeros((nbands, nbands, nk, nk), dtype=complex)

        # Iterate over all possible v, c bands within the specified range
        for iv_index, iv in enumerate(range(min_band, max_band + 1), start=0):
            for ic_index, ic in enumerate(range(min_band, max_band + 1), start=0):
                if iv not in table[:, 1] or ic not in table[:, 2]:
                    continue

                t_v = np.where(table[:, 1] == iv)[0]
                t_c = np.where(table[:, 2] == ic)[0]
                t_vc = [t for t in t_v if t in t_c]

                for it1_subset, it2_subset in product(t_vc, repeat=2):
                    ik = table[it1_subset][0]
                    ip = table[it2_subset][0]
                    # Store in 4D array, adjusting indices for zero-based Python indexing
                    W[iv-1, ic-1, ik - 1, ip - 1] = kernel[it1_subset, it2_subset]

        return W
    
    def get_string(self,mark="="):
        lines = []; app = lines.append
        app( marquee(self.__class__.__name__,mark=mark) )
        app( "kernel mode: RESONANT" )
        app( "number of transitions: %d"%self.ntransitions )
        return '\n'.join(lines)
    
    def get_kernel_indices_bands(self,excitons, bands, iq):
        ''' Given a pair of v and c indices get the t indices of the kernel
        '''
        table  = excitons.table
        nk     = self.lattice.nkpoints
        kernel = self.kernel
        self.consistency_BSE_BSK(excitons)
        
        if bands[0] not in table[:,1] or bands[1] not in table[:,2]:
            raise ValueError('Band indices not matching available transitions')
             
        # Find indices where selected valence band appears
        t_v = np.where(table[:,1]==bands[0])[0]
        # Find indices where selected conduction band appears
        t_c = np.where(table[:,2]==bands[1])[0]
        t_k = np.where(table[:,0]==iq)[0]
        # Among those, find subset of indices where both appear together
        t_vc = [t for t in t_v if (t in t_c and t in t_k) ]

        return t_vc[0]       

    def __str__(self):
        return self.get_string()
