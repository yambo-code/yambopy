# Author: Davide Romanin, FP
#
# This file is part of the yambopy project
#
import os
from yambopy import *
from yambopy.units import *
from itertools import product

class YamboBSEKernelDB(YamboSaveDB):
    """ Read the BSE Kernel database from yambo.
        It reads <t1| K |t2> where K is the kernel and t_i transition indices.
        
        Can only be used if yambo is run with parallel IO.
        
        Only supports "RESONANT" case for BSE calculation.
        TODO: support more cases
    """
    def __init__(self,lattice,excitons,kernel,ker):
        if not isinstance(lattice,YamboLatticeDB):
            raise ValueError('Invalid type for lattice argument. It must be YamboLatticeDB')
        if not isinstance(excitons,YamboExcitonDB):
            raise ValueError('Invalid type for exciton argument. It must be YamboExcitonDB')
        
        self.lattice  = lattice
        self.excitons = excitons
        self.kernel   = kernel
        self.ker      = ker

    @classmethod
    def from_db_file(cls,lattice,excitons,Qpt=1,folder='.'):
        """ initialize this class from a ndb.BS_PAR_Q# file
        """
        filename='ndb.BS_PAR_Q%d'%Qpt
        path_filename = os.path.join(folder,filename)
        if not os.path.isfile(path_filename):
            raise FileNotFoundError("File %s not found in YamboExcitonDB"%path_filename)

        with Dataset(path_filename) as database:
            if 'BSE_RESONANT' in database.variables:
                # Read as transposed since dimension in netCDF are inverted
                reker, imker = database.variables['BSE_RESONANT'][:].T
                ker = reker + imker*I
                # Transform the triangular matrix in a square one
                kernel = np.transpose(np.triu(ker.data))+np.triu(ker.data)
                kernel[np.diag_indices(len(kernel))]*=0.5
            else:
                raise ValueError('Only BSE_RESONANT case supported so far')

        return cls(lattice,excitons,kernel,ker)

    @property
    def ntransitions(self): return len(self.kernel)

    def consistency_BSE_BSK(self):
        """ Check that exciton and kernel dbs are consistent
        """
        if self.excitons.nexcitons != self.ntransitions:
            raise ValueError('Exciton and transition spaces have different dimensions!')
        if self.excitons.ntransitions != self.ntransitions:
            raise ValueError('Mismatch in ntransitions between ExcitonDB and BSEkernelDB!')        

    def get_kernel_exciton_basis(self):
        """ Switch from transition |tq>=|kc,k-qv> to excitonic |lq> basis. 
            In this basis the kernel is diagonal.
            
            <l|K|l> = sum_{t,t'}( <l|t><t|K|t'><t'|l> )
                    = sum_{t,t'}( (A^l_t)^* K_tt' A^l_t' ) 
        
            Here t->kcv according to table from YamboExcitonDB database
        """

        kernel   = self.kernel
        Nstates  = self.ntransitions
        eivs     = self.excitons.eigenvectors
        self.consistency_BSE_BSK()

        # Basis transformation
        kernel_exc_basis = np.zeros(Nstates,dtype=np.complex_)
        for il in range(Nstates):
            kernel_exc_basis[il] = np.dot( np.conj(eivs[il]), np.dot(kernel,eivs[il]) )

        return kernel_exc_basis

    def get_kernel_value_bands(self,bands):#size,bands,yker,yexc):
        """ Get value of kernel matrix elements 
            as a function of k in BZ for fixed c,v bands:
            
                K_cv(k,p) = <ck,vk-q|K|cp,vp-q>
                
            bands = [iv,ic] (NB: enumerated starting from one instead of zero) 
        """
        table  = self.excitons.table
        nk     = self.lattice.nkpoints
        kernel = self.kernel
        self.consistency_BSE_BSK()
        
        if bands[0] not in table[:,1] or bands[1] not in table[:,2]:
            raise ValueError('Band indices not matching available transitions')
             
        # Wcv defined on the full BZ (only a subset will be filled)
        Wcv = np.zeros((nk,nk),dtype=np.complex)
        # Find indices where selected valence band appears
        t_v = np.where(table[:,1]==bands[0])[0]
        # Among those, find subset of indices where conduction band also appears
        t_vc = np.where(table[t_v][:,2]==bands[1])[0]

        # Iterate only on the subset
        for it1_subset, it2_subset in product(t_vc,repeat=2):
            ik = table[it1_subset][0]
            ip = table[it2_subset][0]
            Wcv[ik-1,ip-1] = kernel[it1_subset,it2_subset]

        return Wcv

    def get_string(self,mark="="):
        lines = []; app = lines.append
        app( marquee(self.__class__.__name__,mark=mark) )
        app( "kernel mode: RESONANT" )
        app( "number of transitions: %d"%self.ntransitions )
        return '\n'.join(lines)
    
    def __str__(self):
        return self.get_string()
