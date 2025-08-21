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
    Compute the spin matrix elements ⟨S'|S_z|S⟩ for excitons.

    This function calculates the spin matrix elements for excitons using 
    wavefunctions and spin operators. Easy to use interface. Use 
    compute_exciton_spin() incase you already have those db's

    Parameters
    ----------
    path : str, optional
        Path to the directory containing the `SAVE` folder. Default is `.`.
    bse_dir : str, optional
        Directory containing the BSE (Bethe-Salpeter Equation) data. Default is `'SAVE'`.
    iqpt : int or list/array of ints, optional
        Index or indices of the q-point(s) for which the exciton spin is computed. 
        Default is `1` (Gamma point).
    nstates : int, optional
        Number of exciton states to consider. If `-1`, all states are included. Default is `-1`.
    contribution : {'b', 'e', 'h'}, optional
        Specifies which contribution to compute:
        - `'b'`: Total exciton spin (default).
        - `'e'`: Electron spin only.
        - `'h'`: Hole spin only.
    degen_tol : float, optional
        Degeneracy tolerance for excitons in eV. Default is `1e-2` eV.
    sz : ndarray, optional
        Spin-z operator matrix in the basis of spinor wavefunctions.
        Default is `0.5 * np.array([[1, 0], [0, -1]])`.
    return_dbs_and_spin : bool, optional
        return [latticedb, wfdb, excdb (s), elec_spin_matrix]
        Default is True
    Returns
    -------
    exe_Sz : ndarray
        Spin matrix elements for excitons with shape `(niq, nstates, nstates)`, 
        where `niq` is the number of q-points.
    if return_dbs_and_spin is True, will also return
        [latticedb, wfdb, excdb (s), elec_spin_matrix]
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



