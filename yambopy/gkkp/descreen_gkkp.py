import numpy as np
from yambopy.units import amu2ry

def descreen_el_ph(pseudoscreened_elph,ph_energies,ph_modes,Z,Zval,masses=None):
    """
    This function takes pseudo-screened el-ph matrix elements and "descreens" them returning
    the bare el-ph matrix elements under the approximation that the core static inverse
    dielectric function is proportional to Zval/Z.

    This is an approximate procedure from Eq. (XXX) of REF. (XXX).

    The phonon modes must be the orthonormal phonon modes coming from the diagonalization
    of the dynamical matrix, not the atomic-mass weighted phonon modes which are commonly
    used. In case the latter modes (not normalized) are the only ones available, the atomic
    masses in amu units need to be given as argument.

    Inputs:
    - pseudoscreened_elph: pseudoscreened electron-phonon matrix elements (numpy array)
                           units: same as ph_energies
                           format: [iq][ik][il][is][ib1][ib2]
    - ph_energies: phonon energies (numpy array)
                   units: same as pseudoscreened_elph
                   format: [iq][il]
    - ph_modes: phonon modes (numpy array)
                format: [iq][il][iat][ix]
    - Z: ordered list of atomic numbers (list)
                format: [iat]
    - Zval: ordered list of valence electrons (list)
            format: [iat]
    - masses [optional]: atomic masses (list)
                         units: amu
                         format: [iat]

    NB: ensure that the energy units are consistent between ph_energies and pseudoscreened_elph
        (i.e., hartree, rydberg, etc | QE uses RYDBERG | Yambo uses HARTREE)
    NB: if the masses array is provided, then the eigenvectors will be weighted by the masses
        to restore orthonormality. Remember to give the masses in AMU units.

    Example usage reading the data from LetzElPh output

    ```

    from yambopy import LetzElphElectronPhononDB
    from yambopy.units import ha2ev
    ry2ev = ha2ev/2.              # lelphc uses RYDBERG units
    Z = [42,16,16]                # atomic_numbers (hardcoded here, read with pwin or latticedb)
    Zval = [14,6,6]               # Taken from pseudopotentials, read with qepy utility
    M_amu = [95.95,32.065,32.065] # atomic masses in amu (hardcoded here, read with qepy/yambopy)

    ygpp  = LetzElphElectronPhononDB('PATH/ndb.elph') # dbs generated with "kernel=bare" in lelphc input

    bare_gkkp = descreen_el_ph(ygpp.gkkp,ygpp.ph_energies/ry2ev,ygpp.ph_eigenvectors,Z,Zval,masses=M_amu)

    ```
    """
    na = np.newaxis
    Z = np.array(Z)
    Zval = np.array(Zval)
    if masses is not None:
        masses = amu2ry*np.array(masses)
        ph_modes = ph_modes*np.sqrt( masses[na,na,:,na] )

    # Compute mode overlap matrix weighted by core charge
    L = np.einsum('a,qmax,qnax->qmn',Z/Zval,np.conj(ph_modes),ph_modes,optimize=True)
    # Energy-weighted overlap matrix
    Lw = np.conj(L)*np.sqrt(ph_energies[:,:,na]/ph_energies[:,na,:])
    # Bare matrix elements
    bare_elph =  np.einsum('qmn,qkmsvc->qknsvc',Lw,pseudoscreened_elph,optimize=True)
    
    return bare_elph
