##
## Authors: MN FP
##
import numpy as np
import os
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.bse.exciton_matrix_elements import exciton_X_matelem
from yambopy.bse.rotate_excitonwf import rotate_exc_wf
from tqdm import tqdm
from mpi4py import MPI
comm = MPI.COMM_WORLD
# Variables with rank* are supposed to be different for each process
rank    = comm.Get_rank() # rank id
ntasks = comm.Get_size()

def exciton_phonon_matelem(latdb,elphdb,wfdb,Qrange=[0,1],BSE_dir='bse',BSE_Lin_dir=None,
                           nexc_in=-1,nexc_out=-1,dmat_mode='run',save_files=True,exph_file='Ex-ph.npy',overwrite=False):
    """
    This function calculates the exciton-phonon matrix elements

    - Q is the exciton momentum
    - q is the phonon momentum
    - exc_in represent the "initial" exciton states in the scattering process (at mom. Q)
    - exc_out represents the "final" exciton states in the scattering process (at mom. Q+q)

    MPI support
    -----------
    This is the only function in Yambopy with MPI support. Use it carefully in the following way:
    - read the PARALLEL_EXECUTION.md file
    - make a script that runs this function only (no complicated workflows)
    - the goal is just to have the *.npy exc-ph database printed as a result of the run
    - run as "mpirun -np ntasks python my_input_script_that_calls_this_function.py"
    - use a separate non-MPI script that loads the exc-ph *.npy to proceed with your next calculations

    Parameters
    ----------
    latdb : YamboLatticeDB
        The YamboLatticeDB object which contains the lattice information.
    elphdb : LetzElphElectronPhononDB
        The LetzElphElectronPhononDB object which contains the electron-phonon matrix
        elements.
    wfdb : YamboWFDB
        The YamboWFDB object which contains the wavefunction information.
    BSE_dir : str, optional
        The name of the folder which contains the BSE calculation. Default is 'bse'.
    BSE_Lin_dir : str, optional
        The name of the folder which contains the BSE Q=0 calculation (for optical spectra). Default is BSE_dir.
    Qrange : int list, optional
        Exciton Qpoint index range [iQ_initial, iQ_final] (python counting). Default is [0,0] (Gamma point only).
        Note that the indexing is in full BZ and not in iBZ. See wfc.kBZ to see the kpoints in full BZ
    nexc_in : int, optional
        Number of excitonic states included in Lin calculation. Default is -1 (all).
    nexc_out : int, optional
        Number of excitonic states included in Lout calculation. Default is -1 (all).
    dmat_mode : str, optional
        If 'save', print dmats on .npy file for faster recalculation. If 'load', load from .npy file. Else, calculate Dmats at runtime.
    save_files : bool, optional
        If True, the matrix elements will be saved in .npy file `exph_file`. Default is True.
    overwrite : bool, optional
        If False and `exph_file` is found, the matrix elements will be loaded from file. Default is False.
    """

    # Check if we just need to load
    if os.path.exists(exph_file) and overwrite==False:
        print(f'Loading EXCPH matrices from {exph_file}...')
        exph_mat_loaded = np.load(exph_file)
        return exph_mat_loaded

    # get D matrices
    Dmats = save_or_load_dmat(wfdb,mode=dmat_mode,dmat_file='Dmats.npy')

    # Calculation
    exph_mat = []
    # For each Lin(Q), we call a function that loops over all phonon q
    # and returns <Lin(Q)|dV_ph(q)|Lout(Q+q)> for each q
    for iQ in range(Qrange[0],Qrange[1]):
        Q_in = wfdb.kBZ[iQ]
        #print(f'Calculating EXCPH matrix elements... Q={Q_in}')
        exph_mat.append( exciton_phonon_matelem_iQ(elphdb,wfdb,Dmats,\
                         BSE_dir=BSE_dir,BSE_Lin_dir=BSE_Lin_dir,\
                         Q_in=Q_in,nexc_in=nexc_in,nexc_out=nexc_out) )

    # IO
    if rank==0:
        # single Q-point calculation (suppress axis)
        if len(exph_mat)<2: exph_mat = exph_mat[0]
        # [nQ,nq,nmodes,nexc_in (Qexc),nexc_out (Qexc+q)]
        else:               exph_mat = np.array(exph_mat)
    
        # Save database
        if save_files: 
            if exph_file[-4:]!='.npy': exph_file = exph_file+'.npy'
            print(f'Excph coupling file saved to {exph_file}')
            np.save(exph_file,exph_mat)

    # [WARNING] In MPI execution, rank 0 has the array, all others have None
    return exph_mat

def exciton_phonon_matelem_iQ(elphdb,wfdb,Dmats,BSE_dir,BSE_Lin_dir=None,
                              Q_in=np.zeros(3),nexc_in=-1,nexc_out=-1): 
    """
    This function calculates the exciton-phonon matrix element per Q 

    - Q is the exciton momentum
    - q is the phonon momentum
    - exc_in / Lin represents the "initial" exciton states in the scattering process (at mom. Q)
    - exc_out / Lout represents the "final" exciton states in the scattering process (at mom. Q+q)

    NB: 
        - el-ph couplings are automatically set to convention k->k+q regardless.
        - excitons are instead assumed in convention k-q->k from Yambo.
        - The final results of the exc-ph calculation is INDEPENDENT of conventions.

    Parameters
    ----------
    elphdb : LetzElphElectronPhononDB
        The LetzElphElectronPhononDB object which contains the electron-phonon matrix
        elements.
    wfdb : YamboWFDB
        The YamboWFDB object which contains the wavefunction information.
    BSE_dir : str, optional
        The name of the folder which contains the BSE calculation. Default is 'bse'.
    BSE_Lin_dir : str, optional
        The name of the folder which contains the BSE q=0 calculation (for optical spectra). Default is exdbs[Q].
    Q_in : np.ndarray, optional
        Excitonic momentum in reduced units. Default np.array([0.0,0.0,0.0]) 
    nexc_in : int, optional
        Number of excitonic states included in Lin calculation. Default is -1 (all).
    nexc_out : int, optional
        Number of excitonic states included in Lout calculation. Default is -1 (all).
    """
    latdb = wfdb.ydb
    # Load and rotate Lin(Q)
    if BSE_Lin_dir is None: BSE_Lin_dir = BSE_dir
    Ak = rotate_Akcv_Q(wfdb, Q_in, Dmats, neigs=nexc_in, folder=BSE_Lin_dir)
    if nexc_in == -1: nexc_in = Ak.shape[0]

    # If nexc_out is -1, we need to know it before the MPI chunking
    if nexc_out == -1:
        if rank==0: print("nexc_out is -1, determining the number of excitons from database...")
        idx_BZQ0 = wfdb.kptBZidx(Q_in + elphdb.qpoints[0])
        iQ_iBZ0 = wfdb.ydb.kpoints_indexes[idx_BZQ0]
        filename0 = 'ndb.BS_diago_Q%d' % (iQ_iBZ0+1)
        excdb0 = YamboExcitonDB.from_db_file(wfdb.ydb,filename=filename0,folder=BSE_dir,Load_WF=False, neigs=-1)
        nexc_out = len(excdb0.eigenvalues)
        if rank==0: print(f"nexc_out determined: {nexc_out}")

    # Compute ex-ph
    bse_bnds_range = [wfdb.min_bnd,wfdb.min_bnd + wfdb.nbands]    
    
    #### This loop supports MPI parallelization

    # Basic check
    if ntasks > elphdb.nq:
        if rank==0:
            raise ValueError(f"[MPI][ERROR] You have {ntasks} processes for {elphdb.nq} iterations.")
    # Create "chunks" list with no. of q-iterations per rank
    chunks = [elphdb.nq // ntasks] * ntasks
    for i in range(elphdb.nq % ntasks): chunks[i] += 1 # adjust if not equally div
    # Actual MPI elements (since we are dealing with ndim array, chunks/=counts)
    counts = [ chunk*elphdb.nm*nexc_out*nexc_in for chunk in chunks ] 
    # Indices where each rank must gather in global array
    displs = [sum(counts[:i]) for i in range(ntasks)]
    # Process with rank rank computes rank_chunk iterations
    rank_chunk = chunks[rank]
    # All of the above is done so we can allocate a small local array for each rank
    inside_dims = [elphdb.nm,nexc_out,nexc_in]
    rank_exph = np.zeros([rank_chunk]+inside_dims,dtype=complex)
    # Start of MPI partial loop
    rank_start = sum(chunks[:rank]) # q-index where each rank starts
    # Progress bar just for serial execution
    if ntasks==1: rank_range = tqdm(range(rank_chunk),desc=f"Exc-ph matrix element Q={Q_in}")
    else:         rank_range = range(rank_chunk)
    # Now the calculation starts
    for rank_iq in rank_range:
        iq = rank_start + rank_iq # global q-point index
        # Load el-ph coupling
        ph_eig, elph_mat = elphdb.read_iq(iq,bands_range=bse_bnds_range,convention='standard')
        elph_mat = elph_mat.transpose(1,0,2,4,3)
        # Load and rotate Lout(q+Q)
        Akq = rotate_Akcv_Q(wfdb, Q_in + elphdb.qpoints[iq], Dmats, neigs=nexc_out, folder= BSE_dir)
        # Call the internal generic function to calculate excitonic matrix elements
        rank_exph[rank_iq] = exciton_X_matelem(Q_in, elphdb.qpoints[iq], \
                                              Akq, Ak, elph_mat, wfdb.kBZ, \
                                              contribution='b', diagonal_only=False, \
                                              ktree=wfdb.ktree)

    # Allocate global array on rank 0
    if rank==0: exph_mat = np.zeros([elphdb.nq]+inside_dims,dtype=complex)
    else:       exph_mat = None
    #print(rank,rank_exph.shape)
    #if rank==0: print(counts,displs,exph_mat.shape)
    # Finally, gather call to rank 0
    comm.Gatherv(rank_exph,[exph_mat,counts,displs,MPI.DOUBLE_COMPLEX],root=0)

    if rank==0:
        ## 0.5 for Ry to Ha
        ## change dim order: [nq,nmodes,nexc_in (Qexc),nexc_out (Qexc+q)]
        exph_mat = 0.5 * np.array(exph_mat).transpose(0,1,3,2)
    
    #### End of MPI part 
    return exph_mat

def save_or_load_dmat(wfdb, mode='run', dmat_file='Dmats.npy'):
    """
    Save or load Dmats to/from .npy file `dmat_file` for faster recalculation.

     If mode=='save', print dmats on .npy file for faster recalculation. 
     If mode=='load', load from .npy file. 
     Else, calculate Dmats at runtime.
    """
    if dmat_file[-4:]!='.npy': dmat_file = dmat_file+'.npy'
    if mode=='save':
        if rank==0: print('Saving D matrices...')
        Dmats = wfdb.Dmat()
        if rank==0: np.save(dmat_file,Dmats)
        return Dmats
    elif mode=='load': 
        if rank==0: print('Loading D matrices...')
        if not os.path.exists(dmat_file):
            raise FileNotFoundError(f"Cannot load '{dmat_file}' - file does not exist.")
        Dmats_loaded = np.load(dmat_file)
        return Dmats_loaded
    else:
        return wfdb.Dmat()

def rotate_Akcv_Q(wfdb, Qpt, Dmats, folder, neigs=-1):
    '''
    First load and then rotate exciton coefficients

    wfdb : wavefunction object
    exdbs : previously read list of exciton objects for Lout (used if folder is None)
    Qpt: reduced coordinates in BZ or whatever
    Dmats: precalculated Dmat
    folder: where to load the L exciton states 
    neigs: number of states (used if folder is not None, otherwise it's nexc_out)

    if folder=BSE_Lin_dir and neigs = nexc_in and Qpt = Q --> load Lin
    if folder=BSE_dir and neigs = nexc_out and Qpt = Q+q --> load Lout
    '''
    latdb = wfdb.ydb
    idx_BZQ = wfdb.kptBZidx(Qpt)
    iQ_isymm = latdb.symmetry_indexes[idx_BZQ]
    iQ_iBZ = latdb.kpoints_indexes[idx_BZQ]
    trev  = (iQ_isymm >= len(latdb.sym_car) / (1 + int(np.rint(latdb.time_rev))))
    symm_mat_red = latdb.lat@latdb.sym_car[iQ_isymm]@np.linalg.inv(latdb.lat)
    exe_iQIBZ = wfdb.kpts_iBZ[iQ_iBZ]
    filename = 'ndb.BS_diago_Q%d' % (iQ_iBZ+1)
    
    # Here we load the required exciton data
    excdb = YamboExcitonDB.from_db_file(latdb,filename=filename,folder=folder,\
                                        Load_WF=True, neigs=neigs)
    AQibz = excdb.get_Akcv()
    
    # NM : Add a sanity check to avoid a disastrous consequence
    #      if the user gives wrong bse band indices.
    min_bnd_bse = np.min(excdb.unique_vbands)
    max_bnd_bse = np.max(excdb.unique_cbands)+1
    assert (wfdb.min_bnd == min_bnd_bse) and ((wfdb.min_bnd + wfdb.nbands) == max_bnd_bse), print("[ERROR]: BSE bands mismatch. Given bands range : [%d, %d]. " %(wfdb.min_bnd,wfdb.min_bnd + wfdb.nbands) + "BSE band range found (expected) : [%d %d]" %( min_bnd_bse,max_bnd_bse))
    
    # Here we finally rotate the Akcv coefficients via the internal function
    AQ_rot = rotate_exc_wf(AQibz,symm_mat_red,wfdb.kBZ,exe_iQIBZ,Dmats[iQ_isymm],trev,wfdb.ktree)
    
    return AQ_rot
