import numpy as np
from scipy.special import j0, y0, k0, j1, y1, k1  # Bessel functions from scipy.special
from yambopy.lattice import replicate_red_kmesh, calculate_distances, car_red, modvec
from yambopy.units import alpha,ha2ev, ang2bohr,bohr2ang
import scipy


class CoulombPotentials:
    '''
    class to create Coulomb potential for TB-models. Return values in Hartree
    Example usage
        potentials = CoulombPotentials(ngrid=[...], rlat=[...], tolr=1e-9)
        v2dk_result = potentials.v2dk(kpt1=[...], kpt2=[...], ediel=[...], lc=...)
    '''
    # Constants shared across all methods
    alpha = -(0.0904756) * 10 ** 3
    pi = np.pi

    def __init__(self, ngrid, lattice, v0=0.0, lc=15.0, w=1.0, r0=1.0, tolr=0.001, ediel=[1.0,1.0,1.0]):
        print('''Warning! CoulombPotentials works lat vol in angstrom and kpoints in 1/angstroms.
              ''')
        #lattice is an instance of YamboLatticeDb
        self.ngrid = ngrid
        self.lattice = lattice
        self.rlat = lattice.lat*bohr2ang
        self.tolr = tolr
        self.dir_vol = lattice.lat_vol*bohr2ang**3  #working in angstrom
        self.rec_vol = lattice.rlat_vol*2*np.pi*ang2bohr # reciprocal lattice volume in Bohr
        self.ediel = ediel # ediel(1) Top substrate, ediel(2) \eps_d, ediel(3) Bot substrate     
        self.lc = lc # in angstrom
        self.w = w
        self.r0 = r0
        self.v0 = v0
        self.alpha = -0.0904756*10**3

    def v2dk(self, kpt1, kpt2):
        pass
        #TO-DO implement alat2D
        #constants -> See paper in WantiBexos doc
        alpha1 = 1.76
        alpha2 = 1.0
        alpha3 = 0.0
        lc = self.lc
        ediel = self.ediel
        vsize = np.zeros(3,dtype = np.float64)
        vsize[0] = np.linalg.norm(self.rlat[0])
        vsize[1] = np.linalg.norm(self.rlat[1])
        vsize[2] = np.linalg.norm(self.rlat[2])
        a0 = vsize[0]/2+vsize[1]/2        
        modk = modvec(kpt1, kpt2)
        #compute area of unit cell
        vc = np.linalg.norm(np.cross(self.rlat[0],self.rlat[1]))

        r0 = ((ediel[1] - 1.0) * lc) / (ediel[0] + ediel[2])
        vbz = 1.0 / (np.prod(self.ngrid) * vc)
        gridaux1 = np.float64(self.ngrid[0] * self.ngrid[1])
        auxi = (2.0 * self.pi * r0) / (a0 * np.sqrt(gridaux1))
        ed = (ediel[0] + ediel[2]) / 2.0

        if modk < self.tolr:
            # here factor 2.0*self.pi gets divided by 2*pi
            v2dk = vbz * (self.alpha* a0 * np.sqrt(gridaux1))/(2.*np.pi*ed) * (alpha1 + auxi * alpha2 + alpha3 * auxi**2)
        else:
            v2dk = vbz * (self.alpha/ed) * (1.0 / (modk * (1.0 + r0 * modk)))

        return v2dk*ha2ev

    def vcoul(self, kpt1, kpt2):
        modk = modvec(kpt1, kpt2)
        vbz = self.rec_vol
        ed = self.ediel[1]  # The dielectric constant is set to 1.0
        if modk < self.tolr:
            vcoul = self.v0
        else:
            vcoul = vbz * (self.alpha / ed) * (1.0 / (modk ** 2))

        return vcoul*ha2ev

    def v2dt(self, kpt1, kpt2):
        vc = self.dir_vol
        vbz = 1.0 / (np.prod(self.ngrid) * vc)
        vkpt = np.array(kpt1) - np.array(kpt2)
        gz = abs(vkpt[2])
        gpar = np.sqrt(vkpt[0]**2 + vkpt[1]**2)
        rc = 0.5 * self.rlat[2, 2]
        factor = 1.0  # or 4.0 * self.pi if needed
        modk = modvec(kpt1,kpt2)

        if gpar < self.tolr and gz < self.tolr:
            v2dt = (vbz * self.alpha) * (1/2.0 * rc * rc)
        elif gpar < self.tolr and gz >= self.tolr:
            v2dt = (vbz * self.alpha) * (factor / modk**2) * (1.0 - np.cos(gz * rc) - (gz * rc * np.sin(gz * rc)))
        else:
            aux1 = gz / gpar
            aux2 = gpar * rc
            aux3 = gz * rc
            aux4 = aux1 * np.sin(aux3)
            aux5 = np.cos(aux3)
            v2dt = (vbz * 2*self.alpha) * (factor / modk**2) * (1.0 + (np.exp(-aux2) * (aux4 - aux5)))

        return v2dt*ha2ev

    def v2dt2(self, kpt1, kpt2):
        kpt1_broadcasted = kpt1[:, np.newaxis, :]  # Shape (N1, 1, 3)
        kpt2_broadcasted = kpt2[np.newaxis, :, :]  # Shape (1, N2, 3)
            
        # Compute modk for all kpt1 and kpt2 pairs
        modk = scipy.linalg.norm(kpt1_broadcasted - kpt2_broadcasted, axis=-1)
        
        vbz = 1.0 / (np.prod(self.ngrid) * self.dir_vol)
        lc = self.lc
        vkpt = kpt1_broadcasted - kpt2_broadcasted
            
        # In-plane momentum transfer
        Qxy = np.sqrt(vkpt[..., 0]**2 + vkpt[..., 1]**2)
        # Factor for the potential
        factor = 1.0  # This could also be 4.0 * pi, depending on the model

        # Compute the potential using vectorized operations
        safe_modk = np.where(modk < self.tolr, np.inf, modk)
        v2dt2 = (vbz * self.alpha) * (factor / safe_modk**2) * \
                    (1.0 - np.exp(-0.5 * Qxy * lc) * np.cos(0.5 * lc * vkpt[..., 2]))
        v2dt2 = np.where(modk < self.tolr, 0.0 + 0.j, v2dt2)

        return v2dt2
    

    def v2drk(self, kpt1, kpt2):
        '''
         V(Q) = - e^2 / (2 A_uc \epsilon_0 |Q| F(Q) ) * e^(-w_0 |Q|)
        '''
        w = self.w
        r0 = self.r0
        lc = self.lc
        # Compute the volume of the cell and modulus of the k-point difference
        kpt1_broadcasted = kpt1[:, np.newaxis, :]  # Shape (N1, 1, 3)
        kpt2_broadcasted = kpt2[np.newaxis, :, :]  # Shape (1, N2, 3)
            
        # Compute modk for all kpt1 and kpt2 pairs
        modk = scipy.linalg.norm(kpt1_broadcasted - kpt2_broadcasted, axis=-1)
        
        vbz = 1.0 / (self.rec_vol)

        epar = self.ediel[1]
        et = self.ediel[1]
        eb = self.ediel[0]
        ez = self.ediel[2]
        eta = np.sqrt(epar / ez)
        kappa = np.sqrt(epar * ez)

        pb = (eb - kappa) / (eb + kappa)
        pt = (et - kappa) / (et + kappa)


        aux1 = (1.0 - (pb * pt * np.exp(-2.0 * modk * eta * lc))) * kappa
        aux2 = (1.0 - (pt * np.exp(-eta * modk * lc))) * (1.0 - (pb * np.exp(-eta * modk * lc)))
        aux3 = r0 * modk * np.exp(-modk * w)

        ew = (aux1 / aux2) + aux3       #F(Q)


        safe_modk = np.where(modk < self.tolr, np.inf, modk)
        v2drk = (vbz * self.alpha) * np.exp(-modk * w) * (1.0 / ew) * (1.0 / safe_modk)
        v2drk = np.where(modk < self.tolr, 0.0 + 0.j, v2drk)


        return v2drk*ha2ev
    
    def v0dt(self,kpt1,kpt2):

        vbz = 1.0/(np.prod(self.ngrid)*self.dir_vol)
        modk = modvec(kpt1,kpt2)
        rmin = np.min(np.linalg.norm(self.lattice.lat,axis=1))
        cr = 0.5 * rmin

        factor = 1.0

        if modk < self.tolr:
            v0dt = 1/2.*self.alpha*vbz*cr**2
        else:
            v0dt = self.alpha*vbz*factor/(modk**2)*(1-np.cos(cr*modk))
        return v0dt

    def v2doh(self, kpt1, kpt2):

        w = self.w
        vbz = 1.0/(np.prod(self.ngrid)*self.dir_vol)
        modk = modvec(kpt1,kpt2)
        ez = self.ediel[2]

        if modk < self.tolr:
            v2doh = 0.0
        else :
            v2doh = self.alpha*vbz*np.exp(-w*modk)*(1.0/(ez*modk))
        return v2doh
    
    # def v2dt_vectorized(self, k_diffs):
    #     # this function might be useful to avoid do loops
    #     vc = self.dir_vol
    #     vbz = 1.0 / (np.prod(self.ngrid) * vc)
    #     gz = np.abs(k_diffs[:, 2])
    #     gpar = np.sqrt(k_diffs[:, 0]**2 + k_diffs[:, 1]**2)
    #     rc = 0.5 * self.rlat[2, 2]
    #     factor = 1.0 

    #     modk = self.modvec_vectorized(k_diffs)

    #     v2dt = np.zeros_like(gz)

    #     # Case where gpar and gz are below tolerance
    #     mask_tol = (gpar < self.tolr) & (gz < self.tolr)
    #     v2dt[mask_tol] = (vbz * 2 * self.pi) * (1/8.0 * rc * rc)

    #     # Case where gpar is below tolerance and gz is above tolerance
    #     mask_gpar_tol = (gpar < self.tolr) & (gz >= self.tolr)
    #     v2dt[mask_gpar_tol] = (vbz * 2 * self.pi) * (factor / modk[mask_gpar_tol]**2) * (1.0 - np.cos(gz[mask_gpar_tol] * rc) - (gz[mask_gpar_tol] * rc * np.sin(gz[mask_gpar_tol] * rc)))

    #     # General case
    #     # ~ stands for NOT and | for OR in mask numpy arrays
    #     mask_general = ~(mask_tol | mask_gpar_tol)
    #     aux1 = gz[mask_general] / gpar[mask_general]
    #     aux2 = gpar[mask_general] * rc
    #     aux3 = gz[mask_general] * rc
    #     aux4 = aux1 * np.sin(aux3)
    #     aux5 = np.cos(aux3)
    #     v2dt[mask_general] = (vbz * 2 * self.pi) * (factor / modk[mask_general]**2) * (1.0 + (np.exp(-aux2) * (aux4 - aux5)))

    #     return v2dt * ha2ev

    def modvec_vectorized(self, k_diffs):
        # Vectorized modvec function 
            return np.sqrt(np.sum(k_diffs**2, axis=1))