#
# Authors: MN
#
import os
import numpy as np
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
from yambopy.bse.rotate_excitonwf import rotate_exc_wf
from yambopy.tools.degeneracy_finder import find_degeneracy_evs
from yambopy.symmetries.point_group_ops import get_pg_info, decompose_rep2irrep
from yambopy.symmetries.crystal_symmetries import Crystal_Symmetries
from yambopy.tools.citations import citation

@citation("M. Nalabothula et al. arXiv:2511.21540 (2025)")
def compute_exc_rep(path='.', bse_dir='SAVE', iqpt=1, nstates=-1, degen_tol = 1e-3,
                    degen_rtol=1e-3, symm_tol=1e-2, use_save_symmetries=False):
    """
    Perform a group–theoretical analysis of excitonic wavefunctions
    at a given BSE q-point using the crystal symmetries.

    The routine:
      * Loads lattice, exciton, and electron wave-function databases
      * Identifies positive-energy excitons and degeneracy groups
      * Determines the little group of the chosen q-point using spglib
      * Rotates the exciton wavefunctions under each symmetry
      * Builds the representation matrices of the little group
      * Decomposes each excitonic multiplet into irreducible
        representations of the corresponding point group

    .. note::
       For finite q ≠ 0, use ``Lkind="full"`` (default in Yambo).
       For q = 0, ``Lkind="bar"`` should be used as ``"full"``
       breaks symmetries depending on the gauge.

       For magnetic non-collinear systems, this function falls back to SAVE symmetries.
       This is because, spglib requires magnetic moment  for each atom and there
       is no way we can get it from SAVE as of now. Ofcourse, the price we pay is that
       non-symmorphic symmetries are not taken into account and there by much smaller groups.


    Parameters
    ----------
    path : str, optional
        Working directory containing the ``SAVE`` folder.
        Default is current directory ``'.'``.
    bse_dir : str, optional
        Subdirectory containing the BSE diagonalization database
        (e.g. ``SAVE`` or ``GW_BSE``). Default ``'SAVE'``.
    iqpt : int, optional
        Index of the q-point in the BSE calculation (1-based as in Yambo).
        Default: ``1``.
    nstates : int, optional
        Number of positive-energy exciton states to analyze.
        If ``-1`` (default) all positive states are used.
    degen_tol : float, optional
        Absolute tolerance to detect degeneracies in the exciton energies.
        Default ``1e-3``.
    degen_rtol : float, optional
        Relative tolerance to detect degeneracies in the exciton energies.
        Default ``1e-3``.
    symm_tol : float, optional
        Tolerance to Symmetry operators.
        Default ``1e-2``.
    use_save_symmetries : bool, optional
        If True, uses internal symmetries from SAVE instead of spglib symmetries.

    Outputs
    -------
    No return value.
    Prints to stdout:
        * Identified little group and symmetry operations
        * Degenerate groups of excitons (energies and multiplicities)
        * Irreducible-representation decomposition of each multiplet

    Files Cached in BSE directory (if not present)
    -----------------------------
    ``BS_left_ev_Cache.npy`` :
        Left eigenvectors of the BSE Hamiltonian in case of non-TDA.
    ``Dmat_elec_Cache_spglib.npy`` :
        Electronic representation matrices from spglib symmetry operations.

    Note :
    Currently projective representations are not implemented and left for future.

    Raises
    ------
    AssertionError
        1) If no positive-energy exciton eigenvalues are found.
        2) Projective representations for non-symmorphic symmetires at boundary Q points

    Examples
    --------
    >>> compute_exc_rep(path='..', bse_dir='GW_BSE',
    ...                  iqpt=1, nstates=3, degen_tol=1e-2)

    """
    ## NM : FIX ME, Add a check for the below comment. 
    ## Lkind = "full" for q !=0, as Lkind = "bar" should not be used for finite q
    ## for q = 0, Lking = "bar" is recommaded as "full" will 
    ## break symmetries (depends on the chosen direction
    if nstates == 0 :
        print("Warning : Number of states set to 0. so returning without any output.")
        return None
    # Load the lattice database
    lattice = YamboLatticeDB.from_db_file(os.path.join(path, 'SAVE', 'ns.db1'))
    # NM  : For now fall back to save symmetries as we donot have information on magnetic moments of atoms
    # to get the magnetic symmetries of the crystal.
    if lattice.mag_syms and not use_save_symmetries:
        print("Warning : For magnetic symmetries, the code falls back to symmetries in the SAVE")
        use_save_symmetries = True
    #
    ## Get spglib symmetries to have full symmtries of the crystal even though they are not in yambo base
    # Write basic symmetry information to file
    Rotation_matrices_symm = lattice.sym_car
    frac_trans_symm = np.zeros((len(Rotation_matrices_symm),3),dtype=Rotation_matrices_symm.dtype)
    time_rev = int(np.rint(lattice.time_rev))
    if not use_save_symmetries:
        symm = Crystal_Symmetries(lattice,tol=1e-4)
        symm_info_file_name = "Symmetry_info_excitons.txt"
        print("Writing Basic symmetry information to file : ",symm_info_file_name)
        symm.write_symmetry_info_to_file(symm_info_file_name)
        Rotation_matrices_symm = symm.rotations
        frac_trans_symm = symm.translations
        time_rev = False
    #
    # NM : I am assuming that we can have atmost 100 accidental degeneracies (only for TDA).
    # This is a fine if we are not in continum. In that case, use should increase more states.
    nstates_read = -1
    if nstates >0: nstates_read = nstates + 100
    #
    filename = 'ndb.BS_diago_Q%d' % (iqpt)
    excdb = YamboExcitonDB.from_db_file(lattice, filename=filename,
                                             folder=os.path.join(path, bse_dir),
                                             Load_WF=True, neigs=nstates_read)
    # Load the wavefunction database
    wfdb = YamboWFDB(path=path, latdb=lattice,
                      bands_range=[np.min(excdb.table[:, 1]) - 1,
                    np.max(excdb.table[:, 2])])
    Akcv = excdb.get_Akcv()
    Akcv_left = Akcv
    if Akcv.shape[1] == 2:
        left_ev_path = os.path.join(bse_dir,f'BS_left_ev_Cache_qpt_{iqpt}.npy')
        if os.path.exists(left_ev_path):
            print("Found left eigenvectors, loading ...")
            Akcv_left = np.load(left_ev_path)
        else:
            print("Computing left ev ...")
            Akcv_left = np.linalg.inv(Akcv.reshape(len(Akcv),-1)).conj().T.reshape(Akcv.shape)
            np.save(left_ev_path,Akcv_left)
    #
    eigs = excdb.eigenvalues.real
    sort_idx = np.argsort(eigs)
    eigs = eigs[sort_idx].copy()
    pos_idx = np.where(eigs > 0)
    assert len(pos_idx) >0, "No postive eigenvalues found"
    pos_idx = pos_idx[0][0]
    #
    n_pos = len(eigs[pos_idx:])
    #
    if nstates < 0: nstates = n_pos
    else : nstates = min(nstates, n_pos)
    #
    if nstates < n_pos:
        last = eigs[pos_idx + nstates - 1]
        tol = degen_tol + degen_rtol * abs(last)
        extra = np.sum(np.abs(eigs[pos_idx+nstates:pos_idx+n_pos] - last) < tol)
        if extra: print("Warning: More excitonic states included due to degeneracies.")
        nstates += extra
    #
    Ak_r = Akcv[sort_idx][pos_idx:pos_idx+nstates]
    Ak_l = Akcv_left[sort_idx][pos_idx:pos_idx+nstates].conj()
    exe_ene = eigs[pos_idx:pos_idx+nstates]
    #
    degen_idx = find_degeneracy_evs(exe_ene,atol=degen_tol,rtol=degen_rtol)
    uni_eigs = []
    degen_eigs = []
    for i in degen_idx:
        uni_eigs.append(np.mean(exe_ene[i]))
        degen_eigs.append(len(i))
    uni_eigs = np.array(uni_eigs)
    degen_eigs = np.array(degen_eigs,dtype=int)

    excQpt = excdb.car_qpoint
    # Convert the q-point to crystal coordinates
    excQpt = lattice.lat @ excQpt

    lat_vec = lattice.lat
    lat_vec_inv = np.linalg.inv(lat_vec)
    #
    trev_fac = 1 + int(time_rev)
    # The number of symmetries for which we need representations matrices.
    # We only need representations for spatial symmetries.
    # Note yambo always orders time rev symmeties in second half, even for magnetic systems.
    nsym_spatial = len(Rotation_matrices_symm)//trev_fac
    #
    Dmat_path = os.path.join(bse_dir, 'Dmat_elec_Cache_spglib.npy')
    if use_save_symmetries: Dmat_path = os.path.join(bse_dir, 'Dmat_elec_Cache_SAVE.npy')
    #
    if os.path.exists(Dmat_path):
        print("Dmats found. Loading ....")
        dmats = np.load(Dmat_path)
    else:
        print("Dmats not found. Computing ....")
        dmats = wfdb.Dmat(symm_mat=Rotation_matrices_symm[:nsym_spatial],
                          frac_vec=frac_trans_symm[:nsym_spatial], time_rev=False)
        np.save(Dmat_path,dmats)
    #
    ## print some data about the degeneracies
    print('=' * 40)
    print('Group theory analysis for Q point : (%.6f, %.6f, %.6f)' %(excQpt[0],excQpt[1],excQpt[2]))
    print('*' * 40)

    trace_all_real = []
    trace_all_imag = []
    little_group = []
    #
    #
    for isym in range(nsym_spatial):
        symm_mat = Rotation_matrices_symm[isym]
        symm_mat_red = lat_vec@symm_mat@lat_vec_inv
        #
        Sq_minus_q = np.einsum('ij,j->i', symm_mat_red, excQpt) - excQpt
        #print(Sq_minus_q)
        #diff = Sq_minus_q.copy()
        Sq_minus_q = Sq_minus_q - np.rint(Sq_minus_q)
        ## check if Sq = q
        if np.linalg.norm(Sq_minus_q) > 10**-5: continue
        little_group.append(isym + 1)
        tau_dot_k = np.exp(1j * 2 * np.pi *
                       np.dot(excdb.car_qpoint, frac_trans_symm[isym]))
        #assert(np.linalg.norm(Sq_minus_q)<10**-5)
        rot_Akcv = rotate_exc_wf(Ak_r, symm_mat_red, wfdb.kBZ, excQpt,
                                 dmats[isym], False, ktree=wfdb.ktree)
        rep = tau_dot_k*np.einsum('m...,n...->mn',Ak_l,rot_Akcv,optimize=True)
        # NM : In case of non-trivial projective irrep, we cannot do this,
        # so this fails for boundary Q points, for example (0,0,0.5) in bulk hBN.
        ## Check if this is projective rep and exit as it is yet not implemented.
        G0_tmp = np.einsum('ji,j->i', symm_mat, excdb.car_qpoint) - excdb.car_qpoint
        G0_tau_tmp = G0_tmp.dot(frac_trans_symm[isym])
        G0_tau_tmp = G0_tau_tmp-np.rint(G0_tau_tmp)
        if np.abs(G0_tau_tmp) > 1e-3:
            exit("Error : Projective representations are not implemented")
        #print('Symmetry number : ',isym + 1)
        ## print characters
        irrep_sum = 0
        real_trace = []
        imag_trace = []
        for iirepp in range(len(uni_eigs)):
            idegen = degen_eigs[iirepp]
            idegen2 = irrep_sum + idegen
            trace_tmp = np.trace(rep[irrep_sum:idegen2, irrep_sum:idegen2])
            real_trace.append(trace_tmp.real.round(4))
            imag_trace.append(trace_tmp.imag.round(4))
            irrep_sum = idegen2
        # print('Real : ',real_trace)
        # print('Imag : ',imag_trace)
        trace_all_real.append(real_trace)
        trace_all_imag.append(imag_trace)

    little_group = np.array(little_group, dtype=int)
    #
    pg_label, classes, class_dict, char_tab, irreps = get_pg_info(
        Rotation_matrices_symm[little_group - 1])

    print('Little group : ', pg_label)
    print('Little group symmetries : ', little_group)
    # print class info
    print('Classes (symmetry indices in each class): ')
    req_sym_characters = np.zeros(len(classes), dtype=int)
    class_orders = np.zeros(len(classes), dtype=int)
    for ilab, iclass in class_dict.items():
        print("%16s    : " % (classes[ilab]), little_group[np.array(iclass)])
        req_sym_characters[ilab] = min(iclass)
        class_orders[ilab] = len(iclass)
    print()
    trace_all_real = np.array(trace_all_real)
    trace_all_imag = np.array(trace_all_imag)
    trace = trace_all_real + 1j * trace_all_imag
    trace_req = trace[req_sym_characters, :].T
    print("====== Exciton representations ======")
    print("Energy (eV),  degeneracy  : representation")
    print('-' * 40)
    for i in range(len(trace_req)):
        rep_str_tmp = decompose_rep2irrep(trace_req[i], char_tab, len(little_group),
                                          class_orders, irreps,tol=symm_tol)
        print('%.4f        %9d  : ' % (uni_eigs[i], degen_eigs[i]), rep_str_tmp)
    print('*' * 40)



if __name__ == "__main__":
    compute_exc_rep(path='..', bse_dir='GW_BSE', iqpt=1, nstates=3, degen_tol = 1e-2)

