import os
import numpy as np
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from .exciton_matrix_elements import exciton_X_matelem

def compute_exciton_spin(path='.', bse_dir='SAVE', iqpt=1, nstates=-1, contribution='b',
                         sz=0.5 * np.array([[1, 0], [0, -1]]), 
                         lattice=None, wfdb=None, excdb=None):
    """
    Compute the spin matrix elements <S'|S_z|S> for excitons.

    This function calculates the spin matrix elements for excitons using the
    wavefunctions and spin operators. The spin matrix is computed in the basis
    of exciton states, and off-diagonal elements are included. Diagonalization
    of the matrix in degenerate subspaces is required to obtain the spin values.

    Parameters
    ----------
    path : str, optional
        Path to the directory containing the SAVE folder. Default is '.'.
    bse_dir : str, optional
        Directory containing the BSE data. Default is 'SAVE'.
    iqpt : int or list/array of ints, optional
        Index/indices of the q-point(s) for which the exciton spin is computed. Default is 1 (gamma).
    nstates : int, optional
        Number of exciton states to consider. If -1, all states are included. Default is -1.
    sz : array_like, optional
        Spin-z operator matrix in the basis of spinor wavefunctions. Default is 0.5 * [[1, 0], [0, -1]].
    contribution : char, optional
        'b','e', 'h'. If 'b' total spin is computed. 'e'/'h' for only electron/hole spin. 
         Default is 'b' (total spin)

    if latticedb, wfdb, excdb are not None, then the provided db is used and is not read again
    excdb is None or single YamboExcitonDB or list of YamboExcitonDB. If this is given, then iqpt is ignored
    #
    #
    Returns
    -------
    exe_Sz : ndarray
        Spin matrix elements for excitons with shape (niq, nstates, nstates).

    --------------------------------------------------------------------
    Example:
    --------------------------------------------------------------------
    import numpy as np
    from yambopy.bse.exciton_spin import compute_exciton_spin
    #
    #
    # Total spin
    Sz_exe = compute_exciton_spin(bse_dir='GW_BSE',nstates=2)
    #
    # Only Electron spin
    Sz_exe = compute_exciton_spin(bse_dir='GW_BSE',nstates=2,contribution='e')
    #
    # Only Hole spin
    Sz_exe = compute_exciton_spin(bse_dir='GW_BSE',nstates=2,contribution='h')
    #
    #
    # The first two states are degenerate so, we diagonalize in 2x2 sub block
    w = np.linalg.eigvals(Sz_exe)
    print(w) ## spin values of excitons
    #
    #
    """
    ## Check if it single Q or multiple Q's
    if np.isscalar(iqpt): iqpt = [iqpt]
    else : iqpt = list(iqpt)

    # Load the lattice database
    if not lattice:
        lattice = YamboLatticeDB.from_db_file(os.path.join(path, 'SAVE', 'ns.db1'))

    ##
    if not excdb:
        excdb = []
        for iq in iqpt:
            filename = 'ndb.BS_diago_Q%d' % (iq)
            excdb.append(YamboExcitonDB.from_db_file(lattice, filename=filename,
                                        folder=os.path.join(path, bse_dir),
                                        Load_WF=True, neigs=nstates))
    else :
        if type(excdb) != list:
            excdb = [excdb]

    # Load the wavefunction database
    if not wfdb:
        wfdb = YamboWFDB(path=path, latdb=lattice,
                         bands_range=[np.min(excdb[0].table[:, 1]) - 1,
                        np.max(excdb[0].table[:, 2])])
    #
    # Ensure the calculation is valid only for spinor wavefunctions
    assert wfdb.nspinor == 2, "Makes sense only for nspinor = 2"
    #
    # Compute the spin matrix elements in the BZ
    elec_sz = wfdb.get_spin_m_e_BZ(s_z=sz)
    #
    exe_Sz = []
    for ix in range(len(excdb)):
        ## sanity check
        assert np.min(excdb[ix].table[:, 1]) - 1 == wfdb.min_bnd, \
            "wfdb and exciton db are inconsistant (Bands)"
        ## sanity check
        assert np.max(excdb[ix].table[:, 2]) == wfdb.min_bnd + wfdb.nbands, \
            "wfdb and exciton db are inconsistant (Bands)"

        # Convert the exciton database to kcv format
        excdb[ix].convert_to_kcv()
        # Get the exciton q-point in Cartesian coordinates
        excQpt = excdb[ix].car_qpoint

        # Convert the q-point to crystal coordinates
        excQpt = lattice.lat @ excQpt

        # Compute the exciton spin matrix elements <S'|S_z|S>
        exe_Sz.append(exciton_X_matelem(excQpt, np.array([0, 0, 0]), excdb[ix].eigenvectors,
                               excdb[ix].eigenvectors, elec_sz[None, ...], wfdb.kBZ,
                               diagonal_only=False,contribution=contribution))

    # Print a note about the spin matrix
    #print("Note: This is a spin matrix. Diagonalize the matrix in degenerate subspace to get spin values.")
    
    return np.array(exe_Sz)

