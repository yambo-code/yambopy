import numpy as np
from yambopy.wannier.wann_Gfuncs import GreensFunctions
from yambopy.wannier.wann_io import RMN
from yambopy.wannier.wann_utils import *


class TB_dipoles():
    '''dipoles = 1/(\DeltaE+ieta)*<c,k|P_\alpha|v,k>'''
    def __init__(self , nc, nv, bse_nc, bse_nv, nkpoints, eigv, eigvec, \
                 eta, hlm, T_table, BSE_table, h2peigvec,eigv_diff_ttp, eigvecc_t,eigvecv_t,\
                 h2peigv_vck = None, h2peigvec_vck = None, method = 'real',\
                 rmn = None,ktype='IP'):
        # hk, hlm are TBMODEL hamiltonians
        self.ntransitions = nc*nv*nkpoints
        self.nbsetransitions = bse_nc*bse_nv*nkpoints
        self.nc = nc
        self.nv = nv
        self.nkpoints = nkpoints
        self.eigv = eigv
        self.eigvec = eigvec
        self.nb = nc+nv
        self.bse_nv = bse_nv
        self.bse_nc = bse_nc
        self.bse_nb = bse_nv + bse_nc
        self.offset_nv = nv-bse_nv
        # self.eigvv = eigvv
        # self.eigvc = eigvc
        # self.eigvecv = eigvecv
        # self.eigvecc = eigvecc
        self.eta = eta
        self.hlm = hlm
        self.ktype=ktype
        self.eigv_diff_ttp = eigv_diff_ttp
        self.eigvecc_t = eigvecc_t
        self.eigvecv_t = eigvecv_t
        if self.hlm[0][0][0][0] ==0:
            raise ValueError
        else:
            print(f"hlm not zero: {self.hlm[0][0][0][0]}")
        self.method = method
        if(rmn is not None):
            self.rmn = rmn
            self.method = 'position'
        #T_table = [transition, ik, iv, ic] 
        self.T_table = T_table
        #[nkpoints,3,nbands,nbands]
        # self.dipoles = self._get_dipoles(method)
        # self.d_knm = self._get_dipoles_nm(method)
        if (h2peigvec_vck is not None):
            self.h2peigvec_vck = h2peigvec_vck
            self.h2peigv_vck = h2peigv_vck
            self.h2peigvec = h2peigvec
            # self.dipoles_bse = self._get_dipoles_bse(method)
            self.BSE_table = BSE_table
            if(self.ktype=='IP'):
                print('Running IP dipoles')
                #self._get_dipoles_IP(method=method)
            else:
                print('Running BSE dipoles')
                self._get_dipoles_bse(method=method)
        else:
            print('here')
            #self._get_dipoles(method=method)
        if(self.ktype=='IP'):
            print('Running IP oscillator strength')
            #self._get_osc_strength_IP(method)
        else:
            #self._get_osc_strength_IP(method)
            self._get_osc_strength(method)

    # full dipoles matrix, not only cv, needs adaptation
    def _get_dipoles_nm(self, method):
        if (method == 'real'):
            dipoles = np.zeros((self.nkpoints, self.nb,self.nb,3),dtype=np.complex128)
            for n in range(0, self.nb):
                for m in range(0,self.nb):
                    for ik in range(0,self.nkpoints):
                        # E = self.eigv[ik, n]-self.eigv[ik, m]
                        # GR = GreensFunctions(E,0,self.eta).GR
                        #GA = GreensFunctions(E,0,self.eta).GA
                        # hlm is Bohr*eV
                        dipoles[ik, n, m,0] = np.vdot(self.eigvec[ik,:,n],np.dot(self.hlm[ik,:,:,0],self.eigvec[ik,:,m]))
                        dipoles[ik, n, m,1] = np.vdot(self.eigvec[ik,:,n],np.dot(self.hlm[ik,:,:,1],self.eigvec[ik,:,m]))
                        dipoles[ik, n, m,2] = np.vdot(self.eigvec[ik,:,n],np.dot(self.hlm[ik,:,:,2],self.eigvec[ik,:,m]))

        return dipoles#/(HA2EV**3)
    
    def _get_dipoles(self, method):
        if method == 'real':  # Parallelize over kpoints
            import time
            print("Starting dipole matrix formation.\n")
            t0 = time.time()
            dipoles_kcv = np.zeros((self.nkpoints,self.nb,self.nb,3),dtype=np.complex128)
            for t in range(0,self.ntransitions):
                ik = self.T_table[t][0]
                iv = self.T_table[t][1]
                ic = self.T_table[t][2]
                w = self.eigv[ik,ic]
                E = self.eigv[ik,ic]
                GR = GreensFunctions(w=w, E=E, eta=self.eta).GR
                dipoles_kcv[ik, ic, iv,0] = GR*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,0],self.eigvec[ik,:,iv]))
                dipoles_kcv[ik, ic, iv,1] = GR*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,1],self.eigvec[ik,:,iv]))
                dipoles_kcv[ik, ic, iv,2] = GR*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,2],self.eigvec[ik,:,iv]))

            self.dipoles_kcv = dipoles_kcv
            # # Determine the dimension of hlm
            # dim_hlm = 3 #if np.count_nonzero(self.hlm[:,:,:,2]) > 0 else 2

            # # Extract k, v, c from T_table
            # k_indices, v_indices, c_indices = self.T_table.T

            # # Compute Green's function for all transitions
            # w = self.eigv[k_indices, c_indices]
            # E = self.eigv[k_indices, v_indices]
            # GR = GreensFunctions(w=w, E=E, eta=self.eta).GR

            # # Initialize dipoles array
            # dipoles = np.zeros((self.ntransitions, dim_hlm), dtype=np.complex128)

            # # Prepare eigenvectors
            # eigvec_c = self.eigvec[k_indices, :, c_indices]
            # eigvec_v = self.eigvec[k_indices, :, v_indices]
            # print('shape eigvec', eigvec_c.shape)
            # # Compute dipoles
            # for dim in range(dim_hlm):
            #     # Compute the dot product
            #     print('shape hlm',self.hlm[k_indices,:,:,dim].shape)
            #     dot_product = np.einsum('ijk,ij->ik', self.hlm[k_indices, :, :, dim],eigvec_v)
                
            #     # Compute the vdot and multiply with GR
            #     dipoles[:, dim] = GR * np.einsum('ij,jk->i',np.conjugate(eigvec_c),dot_product)
            #     #* np.sum(np.conjugate(eigvec_c) * dot_product, axis=1)

            # # Reshape dipoles to match your original shape
            # final_dipoles = np.zeros((self.ntransitions, 3), dtype=np.complex128)
            # final_dipoles[np.arange(self.ntransitions), :dim_hlm] = dipoles
            # #self.dipoles_kcv = final_dipoles / (HA2EV ** 3)
            # print("Dipoles matrix computed successfully in serial mode.")
            # print(f"Time for Dipoles matrix formation: {time.time() - t0:.2f}")
        if (method == 'yambo'):
            dipoles = np.zeros((self.ntransitions, self.nkpoints, self.nb,self.nb,3),dtype=np.complex128)
            for n in range(0, self.ntransitions):
                for t in self.T_table:
                    ik = t[0]
                    iv = t[1]
                    ic = t[2]
                    # here I want 1/(E_cv-E_vk) so w=\DeltaE and E = 0 in the call to GFs
                    # E = self.eigv[ik, ic]-self.eigv[ik, iv]
                    GR = GreensFunctions(w=self.eigv[ik, ic], E=self.eigv[ik, iv], eta=self.eta).GR #w - E
                    GA = GreensFunctions(w=self.eigv[ik, ic], E=self.eigv[ik, iv], eta=self.eta).GA #w - E
                    dipoles[n, ik, ic, iv, 0] = (GR+GA)*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,0],self.eigvec[ik,:,iv]))
                    dipoles[n, ik, ic, iv, 1] = (GR+GA)*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,1],self.eigvec[ik,:,iv]))
                    dipoles[n, ik, ic, iv, 2] = (GR+GA)*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,2],self.eigvec[ik,:,iv]))
            self.dipoles = dipoles/(HA2EV**3)  

        if (method== 'v-gauge'):
            print('Warning! velocity gauge not implemented yet')
        if (method== 'r-gauge'):
            print('Warning! position gauge not implemented yet')
        if (method== 'covariant'):
            print('Warning! covariant approach not implemented yet')
           

    def _get_dipoles_bse(self, method):
        if method == 'real':
            import time
            print("Starting BSE dipole matrix formation\n")
            t0 = time.time()
            dipoles_bse_kcv = np.zeros((self.nbsetransitions, self.nkpoints, self.bse_nc,self.bse_nv,3),dtype=np.complex128)
            dipoles_bse_kcv_conj = np.zeros((self.nbsetransitions, self.nkpoints, self.bse_nc,self.bse_nv,3),dtype=np.complex128)            
            for t in range(0,self.nbsetransitions):
                for tp in range(0,self.nbsetransitions):
                    ik = self.BSE_table[tp][0]
                    iv = self.BSE_table[tp][1]
                    ic = self.BSE_table[tp][2]
                    w = self.eigv[ik,ic]
                    E = self.eigv[ik,iv]
                    GR = GreensFunctions(w=w, E=E, eta=self.eta).GR           
                    GA = GreensFunctions(w=w, E=E, eta=self.eta).GA 
                    #dipoles_bse_kck lives in the BSE kernel subset that's why we use the indices ic-self.nv and self.bse_nv-self.nv+iv
                    dipoles_bse_kcv[t, ik, ic-self.nv, iv-self.offset_nv,0] = GR*self.h2peigvec_vck[t,iv-self.offset_nv,ic-self.nv,ik]* \
                        np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,0],self.eigvec[ik,:,iv]))
                    dipoles_bse_kcv[t,ik, ic-self.nv, iv-self.offset_nv,1] = GR*self.h2peigvec_vck[t,iv-self.offset_nv,ic-self.nv,ik]* \
                        np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,1],self.eigvec[ik,:,iv]))
                    dipoles_bse_kcv[t,ik, ic-self.nv, iv-self.offset_nv,2] = GR*self.h2peigvec_vck[t,iv-self.offset_nv,ic-self.nv,ik]* \
                        np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,2],self.eigvec[ik,:,iv]))       
                    dipoles_bse_kcv_conj[t, ik, ic-self.nv, iv-self.offset_nv,0] = GA*self.h2peigvec_vck[t,iv-self.offset_nv,ic-self.nv,ik].conj()* \
                        np.vdot(self.eigvec[ik,:,iv],np.dot(self.hlm[ik,:,:,0],self.eigvec[ik,:,ic]))
                    dipoles_bse_kcv_conj[t,ik, ic-self.nv, iv-self.offset_nv,1] = GA*self.h2peigvec_vck[t,iv-self.offset_nv,ic-self.nv,ik].conj()* \
                        np.vdot(self.eigvec[ik,:,iv],np.dot(self.hlm[ik,:,:,1],self.eigvec[ik,:,ic]))
                    dipoles_bse_kcv_conj[t,ik, ic-self.nv, iv-self.offset_nv,2] = GA*self.h2peigvec_vck[t,iv-self.offset_nv,ic-self.nv,ik].conj()* \
                        np.vdot(self.eigvec[ik,:,iv],np.dot(self.hlm[ik,:,:,2],self.eigvec[ik,:,ic]))                                
            # Determine the dimension of hlm
            #dim_hlm = 3 #if np.count_nonzero(self.hlm[:,:,:,2]) > 0 else 2

            # # Extract k, v, c from T_table
            # k_indices, v_indices, c_indices = self.T_table.T
            # # k_indices = np.array(k_indices)[:, np.newaxis]


            # # Compute Green's function for all transitions
            # E = self.eigv[k_indices, c_indices] - self.eigv[k_indices, v_indices]
            # GR = GreensFunctions(w=E, E=0, eta=self.eta).GR

            # # Initialize dipoles array
            # # dipoles = np.zeros((self.nbsetransitions, dim_hlm), dtype=np.complex128)
            # print(self.nbsetransitions)
            # print(self.bse_nb)
            # print(self.nkpoints)

            # dipoles = np.zeros((self.nbsetransitions, dim_hlm), dtype=np.complex128)
            # # Prepare eigenvectors
            # eigvec_c = self.eigvec[k_indices, :, c_indices]
            # eigvec_v = self.eigvec[k_indices, :, v_indices]

            # # Compute dipoles
            # for dim in range(dim_hlm):
            #     # Compute the dot product
            #     dot_product = np.einsum('ijk,ij->ik', self.hlm[k_indices, :, :, dim],eigvec_v)
                
            #     # Compute the vdot and multiply with GR and h2peigvec
            #     dipoles_cvk[:, :, :, :, dim] = GR * self.h2peigvec[:,:, self.bse_nv-self.nv+v_indices, c_indices-self.nv, k_indices] * np.sum(np.conjugate(eigvec_c) * dot_product, axis=1)

            # Reshape dipoles to match your original shape
            # final_dipoles = np.zeros((self.nbsetransitions, self.bse_nc, self.bse_nv, self.nkpoints, dim_hlm), dtype=np.complex128)
            # final_dipoles[np.arange(self.nbsetransitions), c_indices-self.nv, self.bse_nv-self.nv+v_indices, self.nkpoints,:dim_hlm] = dipoles_cvk

            self.dipoles_bse_kcv = dipoles_bse_kcv
            self.dipoles_bse_kcv_conj = dipoles_bse_kcv_conj

        print("BSE Dipoles matrix computed successfully.")
        print(f"Time for BSE Dipoles matrix formation: {time.time() - t0:.2f}")
        if (method == 'yambo'):
            dipoles = np.zeros((self.ntransitions, self.nkpoints, self.nb,self.nb,3),dtype=np.complex128)
            for n in range(0, self.ntransitions):
                for t in self.T_table:
                    ik = t[0]
                    iv = t[1]
                    ic = t[2]
                    # here I want 1/(E_cv-E_vk) so w=\DeltaE and E = 0 in the call to GFs
                    E = self.eigv[ik, ic]-self.eigv[ik, iv]
                    GR = GreensFunctions(E,0,self.eta).GR
                    GA = GreensFunctions(E,0,self.eta).GA
                    dipoles[n, ik, ic, iv,0] = (GR+GA)*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,0],self.eigvec[ik,:,iv]))
                    dipoles[n, ik, ic, iv,1] = (GR+GA)*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,1],self.eigvec[ik,:,iv]))
                    dipoles[n, ik, ic, iv,2] = (GR+GA)*np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,2],self.eigvec[ik,:,iv]))      
        if (method== 'v-gauge'):
            print('Warning! velocity gauge not implemented yet')
        if (method== 'r-gauge'):
            print('Warning! position gauge not implemented yet')
        if (method== 'covariant'):
            print('Warning! covariant approach not implemented yet')
            #self.dipoles = dipoles Perhaps if method is yambo we want to store dipoles differently                       

    def _get_dipoles_IP(self, method):
        if method == 'real':
            import time
            print("Starting BSE dipole matrix formation.\n")
            t0 = time.time()
            dipoles_bse_kcv = np.zeros((self.nkpoints, self.bse_nc,self.bse_nv,3),dtype=np.complex128)
            for tp in range(0,self.nbsetransitions): 
                ik = self.BSE_table[tp][0]
                iv = self.BSE_table[tp][1]
                ic = self.BSE_table[tp][2]
                w = self.eigv[ik,ic]
                E = self.eigv[ik,iv]
                GR = GreensFunctions(w=w, E=E, eta=self.eta).GR           
                GA = GreensFunctions(w=w, E=E, eta=self.eta).GA 
                #dipoles_bse_kck lives in the BSE kernel subset that's why we use the indices ic-self.nv and self.bse_nv-self.nv+iv
                dipoles_bse_kcv[ik, ic-self.nv, iv-self.offset_nv,0] = GR * \
                    np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,0],self.eigvec[ik,:,iv]))
                dipoles_bse_kcv[ik, ic-self.nv, iv-self.offset_nv,1] = GR * \
                    np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,1],self.eigvec[ik,:,iv]))
                dipoles_bse_kcv[ik, ic-self.nv, iv-self.offset_nv,2] = GR * \
                    np.vdot(self.eigvec[ik,:,ic],np.dot(self.hlm[ik,:,:,2],self.eigvec[ik,:,iv]))            
        self.dipoles_bse_kcv = dipoles_bse_kcv

        print("BSE Dipoles matrix computed successfully.")
        print(f"Time for BSE Dipoles matrix formation: {time.time() - t0:.2f}")
                
    def _get_osc_strength(self,method):
        '''computes osc strength from dipoles
        F_{\alpha, \beta}^{n, BSE} = ( \Sum_c,v,k = \dfrac{A^n_{c,v,k,0} < c,k | P_{\alpha} |v,k >}{E_{c,k} - E{v,k} + i \eta_1} ) * c.c
        '''
        print('Computing oscillator strenght')
        import time
        t0 = time.time()
        dipoles_bse_kcv = self.dipoles_bse_kcv
        dipoles_bse_kcv_conj = self.dipoles_bse_kcv_conj
        F_kcv = np.zeros((self.nbsetransitions, 3, 3), dtype=np.complex128)   



        if (method == 'real'):
            for t in range(0,self.nbsetransitions):
                tmp_F_left = np.zeros((self.nbsetransitions,3), dtype=np.complex128)
                tmp_F_right = np.zeros((self.nbsetransitions,3), dtype=np.complex128)
                for idip in range(0,self.nbsetransitions):
                    ik = self.BSE_table[idip][0]
                    iv = self.BSE_table[idip][1]
                    ic = self.BSE_table[idip][2]
                    factorLx = dipoles_bse_kcv[t, ik, ic-self.nv, iv-self.offset_nv, 0]
                    factorRx = dipoles_bse_kcv_conj[t, ik, ic-self.nv, iv-self.offset_nv, 0] 
                    factorLy = dipoles_bse_kcv[t, ik, ic-self.nv, iv-self.offset_nv, 1]
                    factorRy = dipoles_bse_kcv_conj[t, ik, ic-self.nv, iv-self.offset_nv, 1]
                    factorLz = dipoles_bse_kcv[t, ik, ic-self.nv, iv-self.offset_nv, 2]
                    factorRz = dipoles_bse_kcv_conj[t, ik, ic-self.nv, iv-self.offset_nv, 2]
                    tmp_F_left[t,0] += factorLx
                    tmp_F_left[t,1] += factorLy
                    tmp_F_left[t,2] += factorLz
                    tmp_F_right[t,0] += factorRx
                    tmp_F_right[t,1] += factorRy
                    tmp_F_right[t,2] += factorRz

                F_kcv[t,0,0] = tmp_F_left[t,0]*tmp_F_right[t,0]
                F_kcv[t,0,1] = tmp_F_left[t,0]*tmp_F_right[t,1]    
                F_kcv[t,0,2] = tmp_F_left[t,0]*tmp_F_right[t,2]    
                F_kcv[t,1,0] = tmp_F_left[t,1]*tmp_F_right[t,0]
                F_kcv[t,1,1] = tmp_F_left[t,1]*tmp_F_right[t,1]    
                F_kcv[t,1,2] = tmp_F_left[t,1]*tmp_F_right[t,2]
                F_kcv[t,2,0] = tmp_F_left[t,2]*tmp_F_right[t,0]
                F_kcv[t,2,1] = tmp_F_left[t,2]*tmp_F_right[t,1]    
                F_kcv[t,2,2] = tmp_F_left[t,2]*tmp_F_right[t,2]                                       

                

        if (method== 'v-gauge'):
            print('Warning! velocity gauge not implemented yet')
        if (method== 'r-gauge'):
            print('Warning! position gauge not implemented yet')
        if (method== 'covariant'):
            print('Warning! covariant approach not implemented yet')
        self.F_kcv = F_kcv        
        print(f"Oscillation strenght computed succesfully in {time.time()-t0:.2f}s")

    def _get_osc_strength_IP(self,method):
        '''computes osc strength from dipoles
        '''
        print('Computing oscillator strenght')
        import time
        t0 = time.time()
        dipoles_bse_kcv = self.dipoles_bse_kcv
        F_kcv = np.zeros((self.nkpoints,self.bse_nc,self.bse_nv,3, 3), dtype=np.complex128)   



        if (method == 'real'):
                tmp_F_left = np.zeros((self.nkpoints,self.bse_nc,self.bse_nv,3), dtype=np.complex128)
                tmp_F_right = np.zeros((self.nkpoints,self.bse_nc,self.bse_nv,3), dtype=np.complex128)
                for idip in range(0,self.nbsetransitions):
                    ik = self.BSE_table[idip][0]
                    iv = self.BSE_table[idip][1]
                    ic = self.BSE_table[idip][2]
                    factorLx = dipoles_bse_kcv[ik, ic-self.nv, iv-self.offset_nv, 0]
                    factorRx = factorLx.conj() 
                    factorLy = dipoles_bse_kcv[ik, ic-self.nv, iv-self.offset_nv, 1]
                    factorRy = factorLy.conj() 
                    factorLz = dipoles_bse_kcv[ik, ic-self.nv, iv-self.offset_nv, 2]
                    factorRz = factorLz.conj() 
                    tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,0] = factorLx
                    tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,1] = factorLy
                    tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,2] = factorLz
                    tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,0] = factorRx
                    tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,1] = factorRy
                    tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,2] = factorRz

                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,0,0] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,0]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,0]
                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,0,1] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,0]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,1]    
                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,0,2] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,0]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,2]    
                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,1,0] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,1]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,0]
                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,1,1] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,1]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,1]    
                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,1,2] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,1]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,2]
                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,2,0] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,2]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,0]
                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,2,1] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,2]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,1]    
                    F_kcv[ik,ic-self.nv,iv-self.offset_nv,2,2] = tmp_F_left[ik,ic-self.nv,iv-self.offset_nv,2]*tmp_F_right[ik,ic-self.nv,iv-self.offset_nv,2]                                             

        if (method== 'v-gauge'):
            print('Warning! velocity gauge not implemented yet')
        if (method== 'r-gauge'):
            print('Warning! position gauge not implemented yet')
        if (method== 'covariant'):
            print('Warning! covariant approach not implemented yet')
        self.F_kcv = F_kcv        
        print(f"Oscillation strenght computed succesfully in {time.time()-t0:.2f}s")