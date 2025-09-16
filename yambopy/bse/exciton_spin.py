import os
import numpy as np
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from .exciton_matrix_elements import exciton_X_matelem
from yambopy.tools.degeneracy_finder import find_degeneracy_evs


def compute_exciton_spin(lattice, excdb, wfdb, elec_sz, contribution='b',diagonal=False):
    """
    Compute the spin matrix elements <S'|S_z|S> for excitons.

    This function calculates the spin matrix elements for excitons using the
    wavefunctions and spin operators. The spin matrix is computed in the basis
    of exciton states, and off-diagonal elements are included. Diagonalization
    of the matrix in degenerate subspaces is required to obtain the spin values.

    Parameters
    ----------
    lattice : latticedb
        Lattice database
    excdb : exciton db
        Exciton Database 
    wfdb : wfc db
        wavefunction Database 
    elec_sz : ndarray
        Electron spin matrix elements (nk, nbnds. nbnds).
    contribution : str, optional
        Specifies which contribution to compute:
        - 'b': Total spin (default).
        - 'e': Electron spin only.
        - 'h': Hole spin only.
    diagonal : bool, optional
        If True, only diagonal spin elements are computed. Default is False.

    Returns
    -------
    exe_Sz : ndarray
        Spin matrix elements for excitons with shape (nstates, nstates).
    """
    #
    # Ensure the calculation is valid only for spinor wavefunctions
    assert wfdb.nspinor == 2, "Makes sense only for nspinor = 2"
    #
    # Sanity check
    assert np.min(excdb.table[:, 1]) - 1 == wfdb.min_bnd, \
            "wfdb and exciton db are inconsistant (Bands)"
    ## sanity check
    assert np.max(excdb.table[:, 2]) == wfdb.min_bnd + wfdb.nbands, \
            "wfdb and exciton db are inconsistant (Bands)"
    #
    assert elec_sz.shape == (wfdb.nkBZ, wfdb.nbands, wfdb.nbands)
    # get Akcv
    Akcv = excdb.get_Akcv()
    #
    # Get the exciton q-point in Cartesian coordinates
    excQpt = excdb.car_qpoint
    #
    # Convert the q-point to crystal coordinates
    excQpt = lattice.lat @ excQpt
    #
    # Compute the exciton spin matrix elements <S'|S_z|S>
    exe_Sz = exciton_X_matelem(excQpt, np.array([0, 0, 0]), Akcv,
                               Akcv, elec_sz[None,:,None,...], wfdb.kBZ,
                               diagonal_only=diagonal,contribution=contribution)
    #
    return exe_Sz[0]





def compute_exc_spin_iqpt(path='.', bse_dir='SAVE', iqpt=1,
                          nstates=-1, contribution='b', degen_tol = 1e-2,
                          sz=0.5 * np.array([[1, 0], [0, -1]]),
                          return_dbs_and_spin=True):
    """
    
    
    Description
    -----------
    Compute expectation value of S_z operator for excitons.

    Parameters
    ----------
    path : str, optional
        Path to the directory containing calculation SAVE and BSE folder.
        Default: '.' (current directory)
    bse_dir : str, optional
        Directory containing BSE calculation data. Default: 'SAVE'
    iqpt : int or array-like, optional
        Q-point index or list of Q-point indices to analyze. Default: 1
        (Fortran indexing)
    nstates : int, optional
        Number of excitonic states to consider. Use -1 for all states. Default: -1
    contribution : str, optional
        Which contribution to compute:
        - 'b': both electron and hole (default)
        - 'e': electron only
        - 'h': hole only
    degen_tol : float, optional
        Tolerance for detecting degenerate states. Default: 1e-2
    sz : ndarray, optional
        S_z operator matrix representation. Default: 0.5 * np.array([[1, 0], [0, -1]])
    return_dbs_and_spin : bool, optional
        If True, returns both spin values and database objects. Default: True

    Returns
    -------
    exe_Sz : ndarray
        Array containing S_z expectation values for excitonic states
    dbs_objects : list, optional
        If return_dbs_and_spin=True, returns [lattice, wfdb, excdb, elec_sz] database objects
    Examples
    --------
    Compute the total spin matrix elements for excitons:

    >>> import numpy as np
    >>> from yambopy.bse.exciton_spin import compute_exc_spin_iqpt
    >>> Sz_exe = compute_exc_spin_iqpt(bse_dir='GW_BSE', nstates=2)
    >>> print(Sz_exe)

    Compute only the electron spin contribution:

    >>> Sz_exe = compute_exc_spin_iqpt(bse_dir='GW_BSE', nstates=2, contribution='e')

    Compute only the hole spin contribution:

    >>> Sz_exe = compute_exc_spin_iqpt(bse_dir='GW_BSE', nstates=2, contribution='h')
    """
    #
    ## Check if it single Q or multiple Q's
    if np.isscalar(iqpt): iqpt = [iqpt]
    else : iqpt = list(iqpt)
    # Load the lattice database
    lattice = YamboLatticeDB.from_db_file(os.path.join(path, 'SAVE', 'ns.db1'))
    ## load exbds
    excdb = []
    for iq in iqpt:
        filename = 'ndb.BS_diago_Q%d' % (iq)
        excdb.append(YamboExcitonDB.from_db_file(lattice, filename=filename,
                                                 folder=os.path.join(path, bse_dir),
                                                 Load_WF=True, neigs=nstates))
    # Load the wavefunction database
    wfdb = YamboWFDB(path=path, latdb=lattice,
                      bands_range=[np.min(excdb[0].table[:, 1]) - 1,
                    np.max(excdb[0].table[:, 2])])
    #
    # Compute the spin matrix elements in the BZ
    elec_sz = wfdb.get_spin_m_e_BZ(s_z=sz)
    #
    exe_Sz = []
    for ixdb in excdb:
        smat = compute_exciton_spin(lattice, ixdb,
                                           wfdb, elec_sz,
                                           contribution=contribution,
                                           diagonal=False)
        smat = get_spinvals(smat, ixdb.eigenvalues, atol=degen_tol)
        ss_tmp = []
        for i in smat: ss_tmp = ss_tmp + list(i)
        exe_Sz.append(ss_tmp)
    #
    exe_Sz = np.array(exe_Sz)
    if return_dbs_and_spin : return exe_Sz,[lattice, wfdb, excdb, elec_sz]
    else : return exe_Sz



def get_spinvals(spin_matrix, eigenvalues, atol=1e-3, rtol=1e-3):
    degen_idx = find_degeneracy_evs(eigenvalues,atol=atol, rtol=rtol)
    spins = []
    for id in degen_idx:
        w = np.linalg.eigvals(spin_matrix[id,:][:,id])
        spins.append(w)
    return spins



