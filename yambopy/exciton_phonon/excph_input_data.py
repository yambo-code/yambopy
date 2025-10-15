##
# Authors: FP
##
import numpy as np
import os
from yambopy import YamboLatticeDB,YamboExcitonDB,LetzElphElectronPhononDB,YamboDipolesDB,YamboWFDB
from yambopy.exciton_phonon.excph_matrix_elements import exciton_phonon_matelem
from yambopy.bse.excitondipoles import exc_dipoles_pol

def exc_ph_get_inputs(lat_path,elph_path,bse_path1,mode='PL',bse_path2=None,wf_path=None,dipoles_path=None,nexc_in='all',nexc_out='all',bands_range=[],phonons_range=[],overwrite=False,dmat_file='Dmats.npy',exph_file='Ex-ph.npy',dip_file='exc_dipoles.npy'):
    """
    This functions creates the necessary inputs for exciton-phonon calculations,
    in the format accepted by the related functions.

    It is intended to facilitate the user and to avoid writing long input scripts,
    if directly loading data from databases is the only thing needed.

    Returns:
    ----------
    * Phonon energies [nq,nmodes] (eV)
    * Exciton energies [nq,nexc_out] (eV)
    * Exciton energies [nexc_in] (eV)
    * Exciton-phonon matrix elements [nq,nmodes,nexc_in,nexc_out] (hartree)
    * if mode=='PL': Exciton dipoles [pol,nexc_in] (bohr)

    Parameters
    ----------
    lat_path : string
        Path to lattice directory (SAVE)
    elph_path : string
        Path to elph directory (LetzElPhC)
    bse_path1 : string
        Path to Lout directory (ndb.BS_diago_Q*)
    mode : string, optional
        * 'PL': prepare input for `excph_luminescence` (default)
        * 'life': prepare input for `excph_lifetimes` (TBD)
    bse_path2 : string, optional
        Path to Lin directory. Default bse_path1.
    wf_path : string, optional
        Path to wfc directory (SAVE). Default lat_path.
    dipoles_path : string, optional
        Path to dipoles directory (ndb.dipoles). Needed for mode=='PL'. Default None.
    nexc_in : int, optional
        Number of Lin excitons. Default all.
    nexc_out : int, optional
        Number of Lout excitons. Default all.
    phonons_range : int list [b_i,b_f], optional
        Number of phonon modes included. Python indexing. Right one is excluded.
    phonons_range : int list [ph_i,ph_f], optional
        Number of phonon modes included. Python indexing. Right one is excluded.
    overwrite : bool, optional
        If True, do not read from existing *.npy files and always perform the calculations
    dmat_file, exph_file, dip_file : strings, optional
        Name of .npy auxiliary output files for dmats, exc-ph couplings, and unprojected dipoles
    """
    if bse_path2 is None: bse_path2 = bse_path1
    if wf_path is None:   wf_path = lat_path
    if mode=='PL' and dipoles_path is None:
        raise ValueError('Please specify `dipoles_path` to ndb.dipoles directory')

    # Load lattice
    lattice = YamboLatticeDB.from_db_file(filename=f'{lat_path}/ns.db1')

    # Load electron-phonon energies
    elph    = LetzElphElectronPhononDB(f'{elph_path}/ndb.elph',read_all=False)
    if len(phonons_range)==0: phonons_range=[0,elph.nm]
    ph_energies = elph.ph_energies[:,phonons_range[0]:phonons_range[1]]

    # Load exciton energies (Lout)
    exc_energies = []
    for iq in range(lattice.ibz_nkpoints):
        filename = 'ndb.BS_diago_Q%d' % (iq+1)
        excdb = YamboExcitonDB.from_db_file(lattice,filename=filename,\
                                            folder=bse_path1,Load_WF=False,\
                                            neigs=nexc_out)
        exc_energies.append( excdb.eigenvalues.real )
    exc_energies = np.array(exc_energies)
    exc_energies = exc_energies[lattice.BZ_to_IBZ_indexes,:]

    # Load exciton energies (Lin)
    if bse_path2 is not None:
        excdb_in = YamboExcitonDB.from_db_file(lattice,filename='ndb.BS_diago_Q1',\
                                               folder=bse_path2,Load_WF=False,\
                                               neigs=nexc_in)
        exc_energies_in = excdb_in.eigenvalues.real
    else:
        exc_energies_in = exc_energies[0]

    # Read wavefunctions (this is not needed if exc-ph are already computed)
    wfcs = YamboWFDB(filename='ns.wf',save=wf_path,latdb=lattice,bands_range=bands_range)

    # Calculate and load exciton-phonon matrix elements
    if os.path.isfile(dmat_file) and not overwrite: dmat_mode='load'
    else: dmat_mode='save'
    excph_couplings = exciton_phonon_matelem(lattice,elph,wfcs,BSE_dir=bse_path1,BSE_Lin_dir=bse_path2,\
                                             neigs=max(nexc_in,nexc_out),overwrite=overwrite,\
                                             dmat_mode=dmat_mode,exph_file=exph_file)
    excph_couplings = excph_couplings[:,phonons_range[0]:phonons_range[1],:nexc_in,:nexc_out]
    
    if mode=='PL':
        # Calculate excitonic dipoles (not projected along field)
        exc_dipoles = exc_dipoles_pol(lat_path,dipoles_path=dipoles_path,bse_path=bse_path2,overwrite=overwrite,dip_file=dip_file)
        exc_dipoles = exc_dipoles[:,:nexc_in]

    # Return data
    if mode=='PL':
        return ph_energies, exc_energies, exc_energies_in, excph_couplings, exc_dipoles 
    else: 
        return ph_energies, exc_energies, exc_energies_in, excph_couplings 
    # if mode=='life': TBD
