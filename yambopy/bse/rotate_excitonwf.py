import numpy as np
from yambopy.kpoints import build_ktree, find_kpt

def rotate_exc_wf(Ak, symm_mat_red, kpoints, exe_qpt, dmats, time_rev, ktree=None):
    """
    Rotate the exciton wavefunction Ak using symmetry operations.

    This function applies a symmetry operation to the exciton wavefunction Ak, 
    which is represented in the basis of electronic states. The rotation is 
    performed using the symmetry matrix in reduced coordinates and the 
    corresponding representation matrices.

    Parameters
    ----------
    Ak : array_like
        Exciton wavefunction coefficients with shape (n_exe_states, nk, nc, nv).
        In case of non-TDA, it is (n_exe_states, 2, nk, nc, nv)
    symm_mat_red : array_like
        Symmetry matrix in reduced coordinates with shape (3, 3).
    kpoints : array_like
        K-points in the full Brillouin zone (crystal coordinates) with shape (nk, 3).
    exe_qpt : array_like
        Momentum of the exciton (q-point) in crystal coordinates with shape (3,).
    dmats : array_like
        Representation matrices for the symmetry operation with shape (nk, Rk_band, k_band).
    time_rev : bool
        If True, apply time-reversal symmetry to the wavefunction.
    ktree : object, optional
        Pre-built k-point tree for efficient k-point searching. If not provided, it will be built.

    Returns
    -------
    rot_Ak : ndarray
        Rotated exciton wavefunction coefficients with the same shape as Ak.
    """
    # Initialize the rotated Ak array
    rot_Ak = np.zeros(Ak.shape, dtype=Ak.dtype)
    # Check TDA
    tda = True
    if len(Ak.shape) != 4 : tda = False

    nk, nc, nv = Ak.shape[-3:]
    # Build a k-point tree if not provided
    if not ktree:
        ktree = build_ktree(kpoints)

    # Compute the indices of Rk and Rk - q
    Rkpts = kpoints @ symm_mat_red.T  # Rotated k-points
    k_minus_q = kpoints - exe_qpt[None, :]  # k - q
    idx_Rk = find_kpt(ktree, Rkpts)  # Indices of rotated k-points
    idx_k_minus_q = find_kpt(ktree, k_minus_q)  # Indices of k - q

    # Extract the conduction and valence parts of the representation matrices
    Dcc = dmats[:, nv:, nv:]  # Conduction band part
    Dvv = dmats[idx_k_minus_q, :nv, :nv].conj()  # Valence band part (conjugated)

    # Apply time-reversal symmetry if required
    Ak_tmp = Ak
    if time_rev:
        Ak_tmp = Ak.conj()

    # Rotate the Ak wavefunction using the representation matrices
    if tda :
        rot_Ak[:, idx_Rk, ...] = (Dcc[None, ...] @ Ak_tmp) @ (Dvv.transpose(0, 2, 1)[None, ...])
    else :
        ## rotate the resonant part
        rot_Ak[:, 0, idx_Rk, ...] = (Dcc[None, ...] @ Ak_tmp[:,0,...]) @ (Dvv.transpose(0, 2, 1)[None, ...])
        # Rotate the anti-resonant part
        Dvv = dmats[:, :nv, :nv]  
        Dcc = dmats[idx_k_minus_q, nv:, nv:].conj()  
        rot_Ak[:, 1, idx_Rk, ...] = (Dcc[None, ...] @ Ak_tmp[:,1,...]) @ (Dvv.transpose(0, 2, 1)[None, ...])
    # done return the rotate wfc
    return rot_Ak

