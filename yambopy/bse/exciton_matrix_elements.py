###
###
# This file contains a genenal functions to compute
# < S | O | S'>, where O is an operator.
# for Ex: is O is dV_scf, then these ,matrix elements are 
# ex-ph matrix elements. if O is S_z (spin operator), we get
# spin matrix elements of excitons
###
import numpy as np
from yambopy.kpoints import build_ktree, find_kpt

def exciton_X_matelem(exe_kvec, O_qvec, Akq, Ak, Omn, kpts, contribution='b', diagonal_only=False, ktree=None):
    """
    Compute the exciton matrix elements in the Tamm-Dancoff approximation.

    This function calculates the matrix elements <S (final), k+q | O(q) | S' (initial), k>,
    discarding the third term (disconnected diagram). The calculation is performed in the
    Tamm-Dancoff approximation.

    Parameters
    ----------
    exe_kvec : array_like
        Exciton k-vector in crystal coordinates.
    O_qvec : array_like
        Momentum transfer vector q in crystal coordinates.
    Akq : array_like
        Wavefunction coefficients for k+q (bra wfc) with shape (n_exe_states, nk, nc, nv).
    Ak : array_like
        Wavefunction coefficients for k (ket wfc) with shape (n_exe_states, nk, nc, nv).
    Omn : array_like
        Matrix elements of the operator O in the basis of electronic states with shape (nlambda, nk, nm, nn).
    kpts : array_like
        K-points used to construct the BSE with shape (nk, 3).
    contribution : str, optional
        Specifies the contribution to include in the calculation:
        - 'e' : Only electronic contribution.
        - 'h' : Only hole contribution.
        - 'b' : Both electron and hole contributions (default).
    diagonal_only : bool, optional
        If True, only the diagonal terms are computed. Default is False.

    ktree : KDtree, optional
        If None, will build internally, else use the user provided
    Returns
    -------
    ex_O_mat : ndarray
        The computed exciton matrix elements with shape (nlambda, n_exe_states) if diagonal_only is True,
        or (nlambda, n_exe_states (final), n_exe_states (initial)) if diagonal_only is False.
    """
    # Number of arbitrary parameters (lambda) in the Omn matrix
    nlambda = Omn.shape[0]
    #
    # Shape of the wavefunction coefficients
    n_exe_states, nk, nc, nv = Akq.shape
    #
    # Ensure that the shapes of Akq and Ak match
    assert Akq.shape == Ak.shape, "Wavefunction coefficient mismatch"
    #
    # Ensure that the contribution parameter is valid
    assert contribution in ['b', 'e', 'h'], "Allowed values for contribution are 'b', 'e', 'h'"
    #
    # Build a k-point tree for efficient k-point searching
    if not ktree : ktree = build_ktree(kpts)
    #
    # Find the indices of k-Q-q and k-q in the k-point tree
    idx_k_minus_Q_minus_q = find_kpt(ktree, kpts - O_qvec[None, :] - exe_kvec[None, :])  # k-Q-q
    idx_k_minus_q = find_kpt(ktree, kpts - O_qvec[None, :])  # k-q
    #
    # Extract the occupied and unoccupied parts of the Omn matrix
    Occ = Omn[:, idx_k_minus_q, nv:, nv:]  # Occupied part
    Ovv = Omn[:, idx_k_minus_Q_minus_q, :nv, :nv]  # conduction part
    #
    # Ensure the arrays are C-contiguous to reduce cache misses
    Ak_electron = np.ascontiguousarray(Ak[:, idx_k_minus_q, ...])
    Akq_conj = Akq.reshape(n_exe_states, -1).conj()
    #
    # Initialize the output matrix
    if diagonal_only:
        ex_O_mat = np.zeros((nlambda, n_exe_states), dtype=Ak.dtype)  # (nlambda, final, initial)
    else:
        ex_O_mat = np.zeros((nlambda, n_exe_states, n_exe_states), dtype=Ak.dtype)  # (nlambda, final, initial)
    #
    # Loop over the arbitrary parameters (lambda)
    for il in range(nlambda):
        # Compute the electron contribution
        if contribution == 'e' or contribution == 'b':
            tmp_wfc = Occ[il][None, :, :, :] @ Ak_electron
        #
        # Compute the hole contribution and subtract from the electron contribution
        if contribution == 'h' or contribution == 'b':
            tmp_h = -Ak @ Ovv[il][None, ...]
            if contribution == 'b':
                tmp_wfc += tmp_h
            else:
                tmp_wfc = tmp_h
        #
        # Reshape the temporary wavefunction coefficients
        tmp_wfc = tmp_wfc.reshape(n_exe_states, -1)
        #
        # Compute the final matrix elements
        if diagonal_only:
            ex_O_mat[il] = np.sum(Akq_conj * tmp_wfc, axis=-1)
        else:
            np.matmul(Akq_conj, tmp_wfc.T, out=ex_O_mat[il])
    #
    # Return the computed exciton matrix elements
    return ex_O_mat

