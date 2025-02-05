import numpy as np

HA2EV  = 27.211396132
BOHR2ANG = 0.52917720859
BOHR2M = 5.29177210903e-11  # Bohr radius in meters
ANG2BOHR = 1./BOHR2ANG
HBAR = 1.054571817e-34# J*s => w = E/hbar
EVTOJ = 1.60218e-19  # Conversion factor from eV to J
C = 1./137.035999084  # Speed of light
HA2J = 4.3597482e-18 # Hartree to Joule
AU2FS =  0.02418884254 # atomic unit of time



def fermi_dirac_T(e, T, fermie):
    # atomic units
    kb = 3.1669147805081354e-06
    fermi = 1.0/(np.exp(e-fermie)/(kb*T))
    return fermi

def fermi_dirac(e, fermie):
    """ Vectorized Fermi-Dirac function. """
    # Create a boolean mask for conditions
    greater_than_fermie = e > fermie
    less_or_equal_fermie = e <= fermie

    # Initialize the result array with zeros (default case when e > fermie)
    result = np.zeros_like(e)

    # Apply conditions
    result[less_or_equal_fermie] = 1

    return result

def sort_eig(eigv,eigvec=None):
    "Sort eigenvaules and eigenvectors, if given, and convert to real numbers"
    # first take only real parts of the eigenvalues
    tmpeigv=np.array(eigv.real,dtype=np.float64)
    # sort energies
    args=tmpeigv.argsort()
    eigv=eigv[args]
    if not (eigvec is None):
        eigvec=eigvec[args]
        return (eigv,eigvec)
    return eigv

def find_kpoint_index(klist, kpoint):
    """
    Find the index of a kpoint in a list of kpoints.
    
    Parameters:
    - klist: A list or a NumPy array of kpoints.
    - kpoint: A single kpoint to find in the list.
    
    Returns:
    - Index of the kpoint in the list, or a message if the kpoint is not found.
    """
    # Convert klist to a NumPy array for efficient comparison
    klist_np = np.array(klist)
    kpoint_np = np.array(kpoint)

    # Find the index where kpoint matches in klist
    # We use np.all and np.where to compare each kpoint
    indices = np.where(np.all(klist_np == kpoint_np, axis=1))[0]

    if indices.size > 0:
        return indices[0]
    else:
        print('k-point not found')
        return None

def ensure_shape(array, shape, dtype=None):
    if array.shape != shape:
        raise ValueError(f"Expected shape {shape}, but got {array.shape}")
    if dtype is not None and array.dtype != dtype:
        raise TypeError(f"Expected dtype {dtype}, but got {array.dtype}")
    return array
class ElectricField():
    '''
        Creates a time dependent electric-field
    '''  
    def __init__(self, E0, delta = 0.0235, type ='gaussian', Edir=[1.0,1.0,0.0], omega_light = 1.0 , phi=0, \
                 t0=0):
        # delta is related to the FWHM: full width half maximum with relationship FWHM = 2np.sqrt(2ln(2))/delta
        self.E0 = E0 #amplitude electric field
        self.omega_light = omega_light*EVTOJ/HBAR*1e-15 # Convert to 1/s omega_light is supposed to be in eV and then in femtoseconds
        self.phi = phi # phase
        self.Edir = Edir
        self.type = type
        self.delta = delta
        self.t0 = t0


    def E_t(self, t):
        if (self.type == 'monochromatic'):
            return np.dot(self.Edir,self.E0*np.cos(self.omega_light * t + self.phi))
        if (self.type == 'delta'):
            if (t == self.t0):
                return np.dot(self.Edir,self.E0)
            else:
                return np.dot(self.Edir,0.0)
        if (self.type == 'gaussian'):
            return np.dot(self.Edir, self.E0 * np.sin(self.omega_light * t) * np.exp(-self.delta**2 * (t - self.t0)**2 / 2))

class ChangeBasis():
    '''
    Handles change of basis from Wannier to Bloch and viceversa.
    Np = Nk (number of supercell = number of k-points)
    '''
    def __init__(self, model):
        self.irpos = model.irpos
        self.nrpos = model.nrpos
        self.nb = model.nb
        self.Uknm = model.Uknm
        self.eigvec = model.eigvec
        self.nk = model.nk
        self.kmpgrid = model.mpgrid
        self.k = self.kmpgrid.k
        self.car_kpoints = self.kmpgrid.car_kpoints

    #works don't know why. Probably becuase I am storing each component of the term withn the sum
    def _bloch_to_wann_factor(self, ik,  m, ire):
        k_dot_r = np.dot(self.mpgrid.k, self.irpos.T)
        phase_Uknm_nk = np.zeros((self.nb, self.nb, self.nk),dtype=np.complex128)
        for n in range(0,self.nb):
            phase_Uknm_nk[:,n,ik] =  np.exp(-1j*2*np.pi*k_dot_r[ik, ire])*self.Uknm[ik,n,m]*self.eigvec[ik,:,n]#
        return phase_Uknm_nk

    def bloch_to_wann(self):
        mRe = np.zeros((self.nb, self.nb, self.nrpos),dtype=np.complex128)
        for m in range(0,self.nb):
            for ip in range(0,self.nrpos):
                for n in range(0,self.nb):
                    for ik,k in enumerate(self.mpgrid.k):
                        mRe[:,m,ip] += self._bloch_to_wann_factor(self, ik, m, ip)[:,n,ik]
        return mRe

    def _wann_to_bloch_factor(self,mRe, ik,  n, ire):
        k_dot_r = np.dot(self.mpgrid.k, self.irpos.T)
        phase_Uknm_mRe = np.zeros((self.nb, self.nb, self.nrpos),dtype=np.complex128)
        for m in range(0,self.nb):
            phase_Uknm_mRe[:,m,ire] =  np.exp(1j*2*np.pi*k_dot_r[ik, ire]) *self.Uknm[ik,m,n]*mRe[:,m,ire]#
        return phase_Uknm_mRe

    def wann_to_bloch(self, mRe):
        nk_vec = np.zeros((self.nb, self.nb, self.nk),dtype=np.complex128)
        for ik, k in enumerate(self.mpgrid.k):
            for n in range(0,self.nb):
                for m in range(0, self.nb):
                    for p in range(0,self.nrpos):
                        nk_vec[:,n,ik] += 1/self.nk*self._wann_to_bloch_factor(self, mRe, ik, n, p)[:,m,p]
        nk_vec = nk_vec.transpose(1,0,2)
        return nk_vec