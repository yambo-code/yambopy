import numpy as np
from yambopy.wannier.wann_utils import ensure_shape
class ChernNumber():
    def __init__(self, h2p = None, h = None):
        
        if (h2p is not None): 
            self.h2p = h2p
            self.NX, self.NY, self.NZ = self.h2p.qmpgrid.grid_shape
            self.spacing_x = np.array([1/self.NX, 0.0, 0.0], dtype=np.float128)
            # Create masks for points lying on the planes      
            self.i_x = self.h2p.kmpgrid.find_closest_kpoint(self.spacing_x)
            self.spacing_y = np.array([0.0, 1/self.NY, 0.0], dtype=np.float128)
            self.i_y = self.h2p.kmpgrid.find_closest_kpoint(self.spacing_y)
            self.spacing_z = np.array([0.0, 0.0, 1/self.NZ], dtype=np.float128)
            self.i_z = self.h2p.kmpgrid.find_closest_kpoint(self.spacing_z)        
            self.qpdqx_grid = self.h2p.kplusq_table     
            self.spacing_xy = np.array([1/self.NX, 1/self.NY, 0.0], dtype=np.float128)
            self.i_xy = self.h2p.kmpgrid.find_closest_kpoint(self.spacing_xy)   
            self.spacing_zx = np.array([1/self.NX, 0.0, 1/self.NZ], dtype=np.float128)
            self.i_zx = self.h2p.kmpgrid.find_closest_kpoint(self.spacing_zx)           
            self.spacing_yz = np.array([0.0, 1/self.NY, 1/self.NZ], dtype=np.float128)
            self.i_yz = self.h2p.kmpgrid.find_closest_kpoint(self.spacing_yz)               
            # Extract points on each plane
            self.qpdx_plane = self.qpdqx_grid[:,self.i_x][:,1]#q_grid[qpdqx_grid[:,i_x][1]]
            self.qx_plane   = self.qpdqx_grid[:,self.i_x][:,0]#q_grid[qpdqx_grid[:,i_x][1]]
            self.qpdy_plane = self.qpdqx_grid[:,self.i_y][:,1]
            self.qy_plane   = self.qpdqx_grid[:,self.i_y][:,0]
            self.qpdz_plane = self.qpdqx_grid[:,self.i_z][:,1]
            self.qz_plane   = self.qpdqx_grid[:,self.i_z][:,0]
            self.qpdxy_plane = self.qpdqx_grid[:,self.i_xy][:,1]
            self.qpdyz_plane = self.qpdqx_grid[:,self.i_yz][:,1]
            self.qpdzx_plane = self.qpdqx_grid[:,self.i_zx][:,1]                 

        if (h is not None): self.h = h

    def chern_exc_exc(self, integrand):
        """
        Compute the flux of $\frac{A* \partial A}{\partial q} through the planes x=0, y=0, and z=0.
        electron reference frame
        Parameters:
            integrand (callable): Function to evaluate the integrand. Takes q_grid as input.
            
        Returns:
            dict: Fluxes through planes {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z}.
            it corresponds to the chern_exc_exc number if integrand are the bse eigenvectors
        """
        integrand = ensure_shape(integrand, (self.h2p.nq_double,self.h2p.dimbse, self.h2p.bse_nv, self.h2p.bse_nc, self.h2p.nk),dtype=np.complex128)

        # Evaluate the integrand at the points on each plane
        integrand_x = integrand[self.qx_plane].conj()*(integrand[self.qpdx_plane] - integrand[self.qx_plane])
        integrand_y = integrand[self.qy_plane].conj()*(integrand[self.qpdy_plane] - integrand[self.qy_plane])
        integrand_z = integrand[self.qz_plane].conj()*(integrand[self.qpdz_plane] - integrand[self.qz_plane])

        flux_x = np.sum(integrand_z+integrand_y,axis=(0,2,3,4)) 
        flux_y = np.sum(integrand_z+integrand_y,axis=(0,2,3,4)) 
        flux_z = np.sum(integrand_x+integrand_y,axis=(0,2,3,4)) 

        return {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z}

    def chern_e_e(self, integrand):
        """
        Compute the flux of $\frac{A* A }{\partial q} through the planes x=0, y=0, and z=0.
        electron reference frame
        Parameters:
            integrand (callable): Function to evaluate the integrand. Takes q_grid as input.
            
        Returns:
            dict: Fluxes through planes {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z}.
            it corresponds to the chern_e_e number if integrand are the bse eigenvectors
        """
        integrand = ensure_shape(integrand, (self.h2p.nq_double,self.h2p.dimbse, self.h2p.bse_nv, self.h2p.bse_nc, self.h2p.nk),dtype=np.complex128)
        #conduction bands tensors
        eigvec_ck = self.h2p.eigvec[self.qx_plane, :, :][:,:,np.unique(self.h2p.BSE_table[:,2])]#[:,np.newaxis,:,:]  # Valence bands
        eigvec_ckdqx = self.h2p.eigvec[self.qpdx_plane, :, :][:,:,np.unique(self.h2p.BSE_table[:,2])]#[np.newaxis,:,:,:]  # Valence bands 
        eigvec_ckdqy = self.h2p.eigvec[self.qpdy_plane, :, :][:,:,np.unique(self.h2p.BSE_table[:,2])]#[np.newaxis,:,:,:]  # Valence bands 
        eigvec_ckdqz = self.h2p.eigvec[self.qpdz_plane, :, :][:,:,np.unique(self.h2p.BSE_table[:,2])]#[np.newaxis,:,:,:]  # Valence bands 
        dotcx = np.einsum('jkl,jkp->jlp',np.conjugate(eigvec_ckdqx-eigvec_ck), eigvec_ck) #l index is conjugated       
        dotcy = np.einsum('jkl,jkp->jlp',np.conjugate(eigvec_ckdqy-eigvec_ck), eigvec_ck) #l index is conjugated       
        dotcz = np.einsum('jkl,jkp->jlp',np.conjugate(eigvec_ckdqz-eigvec_ck), eigvec_ck) #l index is conjugated       

        # Evaluate the integrand at the points on each plane
        integrand_x = np.einsum('abcde,abcie, adi -> b',integrand.conj(), integrand, dotcx)/self.h2p.nq_double
        integrand_y = np.einsum('abcde,abcie, adi -> b',integrand.conj(), integrand, dotcy)/self.h2p.nq_double
        integrand_z = np.einsum('abcde,abcie, adi -> b',integrand.conj(), integrand, dotcz)/self.h2p.nq_double
        flux_x = (integrand_z+integrand_y)/2 # divide by 2 becuase I sum z and y
        flux_y = (integrand_z+integrand_y)/2 
        flux_z = (integrand_x+integrand_y)/2 
        return {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z} 
    
    def chern_h_h(self, integrand):
        """
        Compute the flux of $\frac{A* A }{\partial q} through the planes x=0, y=0, and z=0.
        electron reference frame
        Parameters:
            integrand (callable): Function to evaluate the integrand. Takes q_grid as input.
            
        Returns:
            dict: Fluxes through planes {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z}.
            it corresponds to the chern_h_h number if integrand are the bse eigenvectors
        """
        integrand = ensure_shape(integrand, (self.h2p.nq_double,self.h2p.dimbse, self.h2p.bse_nv, self.h2p.bse_nc, self.h2p.nk),dtype=np.complex128)
        (kmqdq_grid, kmqdq_grid_table) = self.h2p.kmpgrid.get_kmqpdq_grid(self.h2p.qmpgrid)
        # get overlaps valence w_{vv'k-q}
        ikminusq = self.h2p.kminusq_table[:, :, 1]
        ikmqpdqx = kmqdq_grid_table[0][:,:,1] #dx
        ikmqpdqy = kmqdq_grid_table[1][:,:,1] #dy
        ikmqpdqz = kmqdq_grid_table[2][:,:,1] #dz
        eigvec_vkmq = self.h2p.eigvec[ikminusq, :, :][:,:,:,np.unique(self.h2p.BSE_table[:,1])]#[:,np.newaxis,:,:]  # Valence bands
        eigvec_vkmqpdqx = self.h2p.eigvec[ikmqpdqx, :, :][:,:,:,np.unique(self.h2p.BSE_table[:,1])]#[np.newaxis,:,:,:]  # Valence bands 
        eigvec_vkmqpdqy = self.h2p.eigvec[ikmqpdqy, :, :][:,:,:,np.unique(self.h2p.BSE_table[:,1])]#[np.newaxis,:,:,:]  # Valence bands 
        eigvec_vkmqpdqz = self.h2p.eigvec[ikmqpdqz, :, :][:,:,:,np.unique(self.h2p.BSE_table[:,1])]#[np.newaxis,:,:,:]  # Valence bands 
        dotvx = np.einsum('ijkl,ijkp->ijlp',np.conjugate(eigvec_vkmqpdqx-eigvec_vkmq), eigvec_vkmq) #np.conjugate(eigvec_vkmqpdqx-eigvec_vkmq)l index is conjugated       
        dotvy = np.einsum('ijkl,ijkp->ijlp',np.conjugate(eigvec_vkmqpdqy-eigvec_vkmq), eigvec_vkmq) #l index is conjugated       
        dotvz = np.einsum('ijkl,ijkp->ijlp',np.conjugate(eigvec_vkmqpdqz-eigvec_vkmq), eigvec_vkmq) #l index is conjugated       
        integrand_x = np.einsum('abcde,abhde, eahc -> b',integrand.conj(), integrand, dotvx)/self.h2p.nq_double
        integrand_y = np.einsum('abcde,abhde, eahc -> b',integrand.conj(), integrand, dotvy)/self.h2p.nq_double
        integrand_z = np.einsum('abcde,abhde, eahc -> b',integrand.conj(), integrand, dotvz)/self.h2p.nq_double
        flux_x = (integrand_z+integrand_y)/2 
        flux_y = (integrand_z+integrand_y)/2 
        flux_z = (integrand_x+integrand_y)/2 
        return {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z}     

    def chern_exc_e(self, integrand):
        """
        Compute the flux of $\frac{A* \partial A * w_cc'k}{\partial q} through the planes x=0, y=0, and z=0.
        electron reference frame
        Parameters:
            integrand (callable): Function to evaluate the integrand. Takes q_grid as input.
            
        Returns:
            dict: Fluxes through planes {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z}.
            it corresponds to the chern_exc_e number if integrand are the bse eigenvectors
        """
        integrand = ensure_shape(integrand, (self.h2p.nq_double,self.h2p.dimbse, self.h2p.bse_nv, self.h2p.bse_nc, self.h2p.nk),dtype=np.complex128)
        
        eigvec_ck    = self.h2p.eigvec[self.qx_plane, :, :][:,:,np.unique(self.h2p.BSE_table[:,2])]#[:,np.newaxis,:,:]  # Valence bands
        eigvec_ckdqx = self.h2p.eigvec[self.qpdx_plane, :, :][:,:,np.unique(self.h2p.BSE_table[:,2])]#[np.newaxis,:,:,:]  # Valence bands 
        eigvec_ckdqy = self.h2p.eigvec[self.qpdy_plane, :, :][:,:,np.unique(self.h2p.BSE_table[:,2])]#[np.newaxis,:,:,:]  # Valence bands 
        eigvec_ckdqz = self.h2p.eigvec[self.qpdz_plane, :, :][:,:,np.unique(self.h2p.BSE_table[:,2])]#[np.newaxis,:,:,:]  # Valence bands 
        dotcx = np.einsum('jkl,jkp->jlp',np.conjugate(eigvec_ckdqx-eigvec_ck), eigvec_ck) #l index is conjugated       
        dotcy = np.einsum('jkl,jkp->jlp',np.conjugate(eigvec_ckdqy-eigvec_ck), eigvec_ck) #l index is conjugated       
        dotcz = np.einsum('jkl,jkp->jlp',np.conjugate(eigvec_ckdqz-eigvec_ck), eigvec_ck) #l index is conjugated       

        # Evaluate the integrand at the points on each plane
        integrand_x = np.einsum('abcde,abcie, adi -> b', integrand.conj(), integrand[self.qpdx_plane] - integrand[self.qx_plane], dotcx)/self.h2p.nq_double \
                    + np.einsum('abcde,abcie, adi -> b',(integrand[self.qpdx_plane] - integrand[self.qx_plane]).conj(), integrand, dotcx)/self.h2p.nq_double
        integrand_y = np.einsum('abcde,abcie, adi -> b', integrand.conj(), integrand[self.qpdy_plane] - integrand[self.qx_plane], dotcy)/self.h2p.nq_double \
                    + np.einsum('abcde,abcie, adi -> b',(integrand[self.qpdy_plane] - integrand[self.qx_plane]).conj(), integrand, dotcy)/self.h2p.nq_double
        integrand_z = np.einsum('abcde,abcie, adi -> b', integrand.conj(), integrand[self.qpdz_plane] - integrand[self.qx_plane], dotcz)/self.h2p.nq_double \
                    + np.einsum('abcde,abcie, adi -> b',(integrand[self.qpdz_plane] - integrand[self.qx_plane]).conj(), integrand, dotcz)/self.h2p.nq_double
        flux_x = (integrand_z+integrand_y)/2 # divide by 2 becuase I sum z and y
        flux_y = (integrand_z+integrand_y)/2 
        flux_z = (integrand_x+integrand_y)/2 
        return {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z} 

    def chern_exc_h(self, integrand):
        """
        Compute the flux of $\frac{A* \partial A * w_vvc'k}{\partial q} through the planes x=0, y=0, and z=0.
        electron reference frame
        Parameters:
            integrand (callable): Function to evaluate the integrand. Takes q_grid as input.
            
        Returns:
            dict: Fluxes through planes {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z}.
            it corresponds to the chern_exc_h number if integrand are the bse eigenvectors
        """
        integrand = ensure_shape(integrand, (self.h2p.nq_double,self.h2p.dimbse, self.h2p.bse_nv, self.h2p.bse_nc, self.h2p.nk),dtype=np.complex128)
        (kmqdq_grid, kmqdq_grid_table) = self.h2p.kmpgrid.get_kmqpdq_grid(self.h2p.qmpgrid)
        ikminusq = self.h2p.kminusq_table[:, :, 1]
        ikmqpdqx = kmqdq_grid_table[0][:,:,1] #dx
        ikmqpdqy = kmqdq_grid_table[1][:,:,1] #dy
        ikmqpdqz = kmqdq_grid_table[2][:,:,1] #dz

        eigvec_vkmq = self.h2p.eigvec[ikminusq, :, :][:,:,:,np.unique(self.h2p.BSE_table[:,1])]#[:,np.newaxis,:,:]  # Valence bands
        eigvec_vkmqpdqx = self.h2p.eigvec[ikmqpdqx, :, :][:,:,:,np.unique(self.h2p.BSE_table[:,1])]#[np.newaxis,:,:,:]  # Valence bands 
        eigvec_vkmqpdqy = self.h2p.eigvec[ikmqpdqy, :, :][:,:,:,np.unique(self.h2p.BSE_table[:,1])]#[np.newaxis,:,:,:]  # Valence bands 
        eigvec_vkmqpdqz = self.h2p.eigvec[ikmqpdqz, :, :][:,:,:,np.unique(self.h2p.BSE_table[:,1])]#[np.newaxis,:,:,:]  # Valence bands 
        dotvx = np.einsum('ijkl,ijkp->ijlp',np.conjugate(eigvec_vkmqpdqx-eigvec_vkmq), eigvec_vkmq) #np.conjugate(eigvec_vkmqpdqx-eigvec_vkmq)l index is conjugated       
        dotvy = np.einsum('ijkl,ijkp->ijlp',np.conjugate(eigvec_vkmqpdqy-eigvec_vkmq), eigvec_vkmq) #l index is conjugated       
        dotvz = np.einsum('ijkl,ijkp->ijlp',np.conjugate(eigvec_vkmqpdqz-eigvec_vkmq), eigvec_vkmq) #l index is conjugated       

        # Evaluate the integrand at the points on each plane
        integrand_x = np.einsum('abcde,abhde, eahc -> b',integrand.conj(), integrand[self.qpdx_plane] - integrand[self.qx_plane], dotvx)/self.h2p.nq_double \
                    + np.einsum('abcde,abhde, eahc -> b', (integrand[self.qpdx_plane] - integrand[self.qx_plane]).conj(), integrand, dotvx)/self.h2p.nq_double
        integrand_y = np.einsum('abcde,abhde, eahc -> b',integrand.conj(), integrand[self.qpdy_plane] - integrand[self.qx_plane], dotvy)/self.h2p.nq_double \
                    + np.einsum('abcde,abhde, eahc -> b', (integrand[self.qpdy_plane] - integrand[self.qx_plane]).conj(), integrand, dotvy)/self.h2p.nq_double
        integrand_z = np.einsum('abcde,abhde, eahc -> b',integrand.conj(), integrand[self.qpdz_plane] - integrand[self.qx_plane], dotvz)/self.h2p.nq_double \
                    + np.einsum('abcde,abhde, eahc-> b', (integrand[self.qpdz_plane] - integrand[self.qx_plane]).conj(), integrand, dotvz)/self.h2p.nq_double
        flux_x = (integrand_z+integrand_y)/2 # divide by 2 becuase I sum z and y
        flux_y = (integrand_z+integrand_y)/2 
        flux_z = (integrand_x+integrand_y)/2 
        return {'x=0': flux_x, 'y=0': flux_y, 'z=0': flux_z} 

    def curvature_plaquette(self):

        # get overlaps valence w_{vv'k-q}
        NX, NY, NZ = self.h2p.qmpgrid.grid_shape
        plaquettes = {
        'x' : np.zeros((NY)*(NZ),dtype=np.complex128),
        'y' : np.zeros((NZ)*(NX),dtype=np.complex128),
        'z' : np.zeros((NX)*(NY),dtype=np.complex128),
        }       

        plaquettes['x'] = np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdy_plane].conj(),self.h2p.h2peigvec[self.qx_plane])  \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdyz_plane].conj(),self.h2p.h2peigvec[self.qpdy_plane]) \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdz_plane].conj(),self.h2p.h2peigvec[self.qpdyz_plane]) \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qx_plane].conj(),self.h2p.h2peigvec[self.qpdz_plane])  
        plaquettes['y'] = np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdz_plane].conj(),self.h2p.h2peigvec[self.qx_plane])  \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdzx_plane].conj(),self.h2p.h2peigvec[self.qpdz_plane]) \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdx_plane].conj(),self.h2p.h2peigvec[self.qpdzx_plane]) \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qx_plane].conj(),self.h2p.h2peigvec[self.qpdx_plane])   
        plaquettes['z'] = np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdx_plane].conj(),self.h2p.h2peigvec[self.qx_plane]) \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdxy_plane].conj(),self.h2p.h2peigvec[self.qpdx_plane]) \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qpdy_plane].conj(),self.h2p.h2peigvec[self.qpdxy_plane]) \
                        * np.einsum('kts, kts -> ks' , self.h2p.h2peigvec[self.qx_plane].conj(),self.h2p.h2peigvec[self.qpdy_plane])                                      
        return plaquettes