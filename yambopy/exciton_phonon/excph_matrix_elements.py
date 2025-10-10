##
## Authors: MN (FP adapted)
##
import numpy as np
import os
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.bse.exciton_matrix_elements import exciton_X_matelem
from yambopy.bse.rotate_excitonwf import rotate_exc_wf
from tqdm import tqdm

def exciton_phonon_matelem(latdb,elphdb,wfdb,Qrange=[0,0],BSE_dir='bse',BSE_Lin_dir=None,neigs=-1,dmat_mode='run',save_files=True,exph_file='Ex-ph.npy',overwrite=False):
    """
    This function calculates the exciton-phonon matrix elements

    - Q is the exciton momentum
    - q is the phonon momentum
    - exc_in represent the "initial" exciton states in the scattering process (at mom. Q)
    - exc_out represents the "final" exciton states in the scattering process (at mom. Q+q)

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
    neigs : int, optional
        Number of excitonic states included in calculation. Default is -1 (all).
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

    # Load exc dbs
    exdbs = []
    for ik in range(wfdb.nkpoints):
        filename = 'ndb.BS_diago_Q%d' % (ik+1)
        excdb = YamboExcitonDB.from_db_file(latdb,filename=filename,folder=BSE_dir,\
                                            Load_WF=True, neigs=neigs)
        exdbs.append(excdb)

    # get D matrices
    Dmats = save_or_load_dmat(wfdb,mode=dmat_mode,dmat_file='Dmats.npy')

    # Calculation
    print('Calculating EXCPH matrix elements...')
    exph_mat = []
    for iQ in tqdm(range(Qrange[0],Qrange[1]+1)):
        exph_mat.append( exciton_phonon_matelem_iQ(latdb,elphdb,wfdb,exdbs,Dmats,\
                                                   BSE_Lin_dir=BSE_Lin_dir,Qexc=iQ,neigs=neigs) )
    # IO
    if len(exph_mat)<2: exph_mat = exph_mat[0] # single Q-point calculation (suppress axis)
    else:               exph_mat = np.array(excph_mat) #[nQ,nq,nmodes,nexc_in (Qexc),nexc_out (Qexc+q)]
    
    if save_files: 
        if exph_file[-4:]!='.npy': exph_file = exph_file+'.npy'
        print(f'Excph coupling file saved to {exph_file}')
        np.save(exph_file,exph_mat)

def exciton_phonon_matelem_iQ(latdb,elphdb,wfdb,exdbs,Dmats,BSE_Lin_dir=None,Qexc=0,neigs=-1,dmat_mode='run'): 
    """
    This function calculates the exciton-phonon matrix element per Q 

    - Q is the exciton momentum
    - q is the phonon momentum
    - exc_in represent the "initial" exciton states in the scattering process (at mom. Q)
    - exc_out represents the "final" exciton states in the scattering process (at mom. Q+q)

    Parameters
    ----------
    latdb : YamboLatticeDB
        The YamboLatticeDB object which contains the lattice information.
    elphdb : LetzElphElectronPhononDB
        The LetzElphElectronPhononDB object which contains the electron-phonon matrix
        elements.
    wfdb : YamboWFDB
        The YamboWFDB object which contains the wavefunction information.
    exdbs : YamboExcitonDB list
        List of Q+q YamboExcitonDB objects containing the BSE calculation
    BSE_Lin_dir : str, optional
        The name of the folder which contains the BSE q=0 calculation (for optical spectra). Default is exdbs[Q].
    Qexc : int, optional
        Excitonic momentum index (python counting). Default is 0, i.e. Gamma point.
    neigs : int, optional
        Number of excitonic states included in calculation. Default is -1 (all).
    """
    
    # Exc(in) momentum
    Q_in = wfdb.kpts_iBZ[Qexc]

    # Determine Lkind(in)
    if BSE_Lin_dir is None: 
        Ak = exdbs[Qexc].get_Akcv()
    else:
        filename = 'ndb.BS_diago_Q%d' % (Qexc+1)
        excdbin = YamboExcitonDB.from_db_file(latdb,filename=filename,folder=BSE_Lin_dir,\
                                              Load_WF=True, neigs=neigs)
        Ak = excdbin.get_Akcv()
    
    # Compute ex-ph
    exph_mat = []
    for iq in range(elphdb.nq):
        ph_eig, elph_mat = elphdb.read_iq(iq,convention='standard')
        elph_mat = elph_mat.transpose(1,0,2,4,3)
        #
        idx_BZq = wfdb.kptBZidx(elphdb.qpoints[iq])
        iq_isymm = latdb.symmetry_indexes[idx_BZq]
        iq_iBZ = latdb.kpoints_indexes[idx_BZq]
        trev  = (iq_isymm >= len(latdb.sym_car) / (1 + int(np.rint(latdb.time_rev))))
        symm_mat_red = latdb.lat@latdb.sym_car[iq_isymm]@np.linalg.inv(latdb.lat)
        #
        exe_iqvec = wfdb.kpts_iBZ[iq_iBZ]
        #
        Akq = rotate_exc_wf(exdbs[iq_iBZ].get_Akcv(),symm_mat_red,wfdb.kBZ,exe_iqvec,Dmats[iq_isymm],trev,wfdb.ktree)
        #
        tmp_exph = exciton_X_matelem(Q_in, elphdb.qpoints[iq], \
                                     Akq, exdbs[0].get_Akcv(), elph_mat, wfdb.kBZ, \
                                     contribution='b', diagonal_only=False, ktree=wfdb.ktree)
        exph_mat.append(tmp_exph)

    ## 0.5 for Ry to Ha
    exph_mat = 0.5 * np.array(exph_mat).transpose(0,1,3,2) #[nq,nmodes,nexc_in (Qexc),nexc_out (Qexc+q)]

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
        print('Saving D matrices...')
        Dmats = wfdb.Dmat()
        np.save(dmat_file,Dmats)
        return Dmats
    elif mode=='load': 
        print('Loading D matrices...')
        if not os.path.exists(dmat_file):
            raise FileNotFoundError(f"Cannot load '{dmat_file}' - file does not exist.")
        Dmats_loaded = np.load(dmat_file)
        return Dmats_loaded
    else:
        return wfdb.Dmat()

