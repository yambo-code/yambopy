### Compute real space exction wavefunction when hole/electron is fixed.
import numpy as np
from yambopy.kpoints import build_ktree, find_kpt
from tqdm import tqdm
import os

## Usage
"""
import numpy as np
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB
import os

iqpt = 1 # qpt index of exciton

# load lattice db
lattice = YamboLatticeDB.from_db_file(os.path.join('.', 'SAVE', 'ns.db1'))

# load exciton db
# note in case too many excitons, load only first few with `neigs' flag
# DO NOT forget to include all degenerate states when giving neigs flag !
#
filename = 'ndb.BS_diago_Q%d' % (iqpt)
excdb = YamboExcitonDB.from_db_file(lattice, filename=filename,
                                    folder=os.path.join('.', 'GW_BSE'),
                                    neigs = 20)

#Load the wavefunction database
wfdb = YamboWFDB(path='.', latdb=lattice,
                      bands_range=[np.min(excdb.table[:, 1]) - 1,
                      np.max(excdb.table[:, 2])])

## plot the exciton wavefunction with hole fixed at [0,0,0] (crystal coordinates)
# in a [1,1,1] supercell with 80 Ry wf cutoff
# I want to set the degeneracy threshold to 0.01 eV
# For example I want to plot the 3rd exciton, so iexe = 2 (python indexing )
#
excdb.real_wf_to_cube(iexe=2, wfdb=wfdb, fixed_postion=[0,0,0], supercell=[1,1,1],
                        degen_tol=0.01, wfcCutoffRy=80, fix_particle='h')

## fix_particle = 'e' if you want to fix electron and plot hole density
## .cube will be dumped and use vesta to visualize it !
"""

##

def ex_wf2Real(Akcv, Qpt, wfcdb, bse_bnds, fixed_postion,
               fix_particle='h', supercell=[1,1,1], wfcCutoffRy=-1,
               block_size=256):
    """
        Compute real-space exciton wavefunction when hole/electron is fixed.

        This is the main interface function that handles both resonant and anti-resonant parts
        of the exciton wavefunction in the case of non-TDA calculations.

    Args:
        Akcv (numpy.ndarray): Exciton wavefunction coefficients with shape:
                              - (Nstates,k,c,v) for TDA
                              - (Nstates,2,k,c,v) for non-TDA (2 for resonant/anti-resonant)
        Qpt (numpy.ndarray): Q-point of exciton in crystal coordinates
        wfcdb (YamboWFDB): Wavefunction database
        bse_bnds (list): Band range used in BSE [min_band, max_band] (python indexing)
        fixed_postion (list): Position of fixed particle in crystal coordinates
        fix_particle (str): 'e' to fix electron, 'h' to fix hole (default)
        supercell (list): Supercell dimensions [nx,ny,nz]
        wfcCutoffRy (float): Wavefunction cutoff in Rydberg (-1 for no cutoff)
        block_size (int): Block size for memory-efficient computation.
        ## choosing lowe block_size will slight lower the memory requirment but also less faster

    Returns:
        tuple: (supercell_latvecs, atom_nums, atom_pos, exe_wfc_real)
               - supercell_latvecs: Supercell lattice vectors (3,3)
               - atom_nums: Atomic numbers.
               - atom_pos: Atomic positions in cartisian units
               - exe_wfc_real: Real-space exciton wavefunction (nstates, nspinor_electron, nspinor_hole,
                                                            Nx_grid, Ny_grid, Nz_grid)
    """
    if len(Akcv.shape) == 5:
        #nonTDA
        ## first the resonat part
        supercell_latvecs,atom_nums,atom_pos,exe_wfc_real = \
            ex_wf2Real_kernel(Akcv[:,0], Qpt, wfcdb, bse_bnds, fixed_postion,
            fix_particle=fix_particle, supercell=supercell,
            wfcCutoffRy=wfcCutoffRy, block_size=block_size,
                              ares=False, out_res=None)
        # add to antiresonant part
        supercell_latvecs,atom_nums,atom_pos,exe_wfc_real = \
            ex_wf2Real_kernel(Akcv[:,1], Qpt, wfcdb, bse_bnds, fixed_postion,
            fix_particle=fix_particle, supercell=supercell,
            wfcCutoffRy=wfcCutoffRy, block_size=block_size,
                              ares=True, out_res=exe_wfc_real)
    else:
        assert len(Akcv.shape) == 4, "wrong Akcv dimensions"
        supercell_latvecs,atom_nums,atom_pos,exe_wfc_real = \
            ex_wf2Real_kernel(Akcv, Qpt, wfcdb, bse_bnds, fixed_postion,
            fix_particle=fix_particle, supercell=supercell,
            wfcCutoffRy=wfcCutoffRy, block_size=block_size,
                              ares=False, out_res=None)

    return supercell_latvecs,atom_nums,atom_pos,exe_wfc_real



def ex_wf2Real_kernel(Akcv, Qpt, wfcdb, bse_bnds, fixed_postion,
               fix_particle='h', supercell=[1,1,1], wfcCutoffRy=-1,
               block_size=256, ares=False, out_res=None):
    """
        Core kernel function for computing real-space exciton wavefunction.

        Computes either:
        - Akcv * psi_{kc}(r_e) * (psi_{k-Q,v}(r_h))^* (resonant part)
        - Akcv * psi_{kv}(r_e) * (psi_{k-Q,c}(r_h))^* (anti-resonant part)

        Note: For density, one must compute the absolute value squared.

    Args:
        Akcv (numpy.ndarray): Exciton coefficients [nstates,nk,nc,nv]
        Qpt (numpy.ndarray): Q-point in crystal coordinates [3]
        wfcdb (YamboWFDB): Wavefunction database
        bse_bnds (list): BSE band range used in bse [nb1, nb2]. fortran indexing. i.e index starts from 1.
                        i.e nb1  and nb2 are same indices used in yambo input i.e
                        % BSEBands nb1 | nb2 % in yambo input
        fixed_postion (list): Fixed particle position in crystal coords [3]
        fix_particle (str): 'e'=fix electron, 'h'=fix hole (default)
        supercell (list): Supercell dimensions [nx,ny,nz]
        wfcCutoffRy (float): Wavefunction cutoff in Rydberg (-1= full cutoff)
        block_size (int): Memory block size for computation. is a postive integer,
                        the default is 256 which is very good but uses more memory.
    # decrease it when you run into memory issues
        ares (bool): If True, compute anti-resonant part
        out_res (numpy.ndarray): Adds to this array
                                 and is returned instead of internally creating.
                                 Make sure it its consistant. (no internal checking done)

    Returns:
        tuple: (supercell_latvecs, atom_nums, atom_pos, exe_wfc_real)
        Natoms = (natom_in_unit_call * Nsupercell) + 1 (+1 due to hole/electron)
        - supercell_latvecs (numpy.ndarray): Supercell lattice vectors [3,3]
        - atom_nums (numpy.ndarray): Atomic numbers [Natoms]
        - atom_pos (numpy.ndarray): Atomic positions in cartisian units [Natoms,3]
        - exe_wfc_real (numpy.ndarray): Wavefunction in real space
          [nstates, nspinor_electron, nspinor_hole, FFTx, FFTy, FFTz]
    """
    #
    #
    if block_size < 1:
        print('Warning: Wrong block_size. setting to 1')
        block_size = 1
    #
    # Convert them to 
    for i in range(3):
        if supercell[i]%2 == 0:
            print('Warning : Even supercell given, so increasing'
                  + ' supercell size along %d direction by 1'%(i+1))
            supercell[i] = supercell[i] +1
    #
    #
    fix_particle = fix_particle.lower()
    bse_bnds = [min(bse_bnds)-1,max(bse_bnds)]
    assert bse_bnds[0] >= wfcdb.min_bnd, \
        "%d is used in bse but not found in wfcdb, load more wfcs" %(bse_bnds[0]+1)
    assert bse_bnds[1] <= wfcdb.min_bnd + wfcdb.nbands, \
        "%d is used in bse but not found in wfcdb, load more wfcs" %(wfcdb.min_bnd + wfcdb.nbands)
    bse_bnds = np.array(bse_bnds,dtype=int)-wfcdb.min_bnd
    
    nstates, nk, nc, nv = Akcv.shape
    
    kpt_idx = wfcdb.ydb.kpoints_indexes
    sym_idx = wfcdb.ydb.symmetry_indexes
    nkBZ = len(sym_idx)
    assert nc + nv == bse_bnds[1]-bse_bnds[0], "Band mismatch"
    assert nk == nkBZ, "kpoint mismatch"
    #
    
    # ## Bring the positon of hole/electron in the centre unit cell, i.e make in between [0,1)
    # fixed_postion = np.array(fixed_postion) - np.floor(fixed_postion)
    # fixed_postion = (fixed_postion + 1e-6)%1
    # fixed_postion = fixed_postion - 1e-6
    fixed_postion = np.array(fixed_postion)

    fixed_postion += np.array(supercell)//2


    hole_bnds = [bse_bnds[0],bse_bnds[0]+nv]
    elec_bnds = [bse_bnds[0]+nv,bse_bnds[1]]
    #
    if ares:
        # if anti-resonant part, we swap the electron and hole bands
        # first transpose c,v dimensions
        Akcv = Akcv.transpose(0,1,3,2)
        tmp  = hole_bnds
        hole_bnds = elec_bnds
        elec_bnds = tmp
        nstates, nk, nc, nv = Akcv.shape

    lat_vec = wfcdb.ydb.lat.T
    blat = np.linalg.inv(lat_vec)
    gvecs_iBZ_idx = []
    fft_box = np.zeros(3,dtype=int)
    #
    
    for ik in range(len(wfcdb.gvecs)):
        idx_gvecs_tmp = np.arange(wfcdb.ngvecs[ik],dtype=int)
        if wfcCutoffRy > 0:
            tmp_gvecs = np.linalg.norm((wfcdb.gvecs[ik, :wfcdb.ngvecs[ik], :] 
                                        + wfcdb.kpts_iBZ[ik][None,:])@blat,axis=-1)
            idx_tmp = tmp_gvecs < np.sqrt(wfcCutoffRy)
            idx_gvecs_tmp = idx_gvecs_tmp[idx_tmp].copy()
        #
        gvecs_iBZ_idx.append(idx_gvecs_tmp)
        ## Get the fft box 
        min_fft_idx = np.min(wfcdb.gvecs[ik, :wfcdb.ngvecs[ik], :][idx_gvecs_tmp] , axis=0)
        max_fft_idx = np.max(wfcdb.gvecs[ik, :wfcdb.ngvecs[ik], :][idx_gvecs_tmp] , axis=0)
        assert np.min(max_fft_idx) >= 0 and np.max(min_fft_idx) < 0, "Invalid G-vectors"
        for i in range(3):
            fft_box[i] = max([fft_box[i], max_fft_idx[i] - min_fft_idx[i] + 3])
    
    # Compute nstates, nk, Nx, Ny, Nz object
    if out_res is None : print("Wfc FFT Grid : ",fft_box[0], fft_box[1], fft_box[2])
    #
    ktree = build_ktree(wfcdb.kBZ)
    #
    nspinorr = wfcdb.nspinor
    if out_res is not None:
        exe_wfc_real = out_res.reshape(nstates, nspinorr, nspinorr,
                             supercell[0],fft_box[0],
                             supercell[1],fft_box[1],
                             supercell[2],fft_box[2])
    else:
        exe_wfc_real = np.zeros((nstates, nspinorr, nspinorr,
                             supercell[0],fft_box[0],
                             supercell[1],fft_box[1],
                             supercell[2],fft_box[2]),
                            dtype=np.complex64)
    #
    Lx = np.arange(supercell[0],dtype=int)
    Ly = np.arange(supercell[1],dtype=int)
    Lz = np.arange(supercell[2],dtype=int)
    Lsupercells = np.zeros((supercell[0],supercell[1],supercell[2],3),dtype=int)
    Lsupercells[...,0], Lsupercells[...,1], Lsupercells[...,2] = np.meshgrid(Lx,Ly,Lz,indexing='ij')
    #
    FFTxx = np.fft.fftfreq(fft_box[0])
    FFTyy = np.fft.fftfreq(fft_box[1])
    FFTzz = np.fft.fftfreq(fft_box[2])
    FFTxx = FFTxx - np.floor(FFTxx)
    FFTyy = FFTyy - np.floor(FFTyy)
    FFTzz = FFTzz - np.floor(FFTzz)
    FFFboxs = np.zeros((fft_box[0],fft_box[1],fft_box[2],3))
    FFFboxs[...,0], FFFboxs[...,1], FFFboxs[...,2] = np.meshgrid(FFTxx,FFTyy,FFTzz,indexing='ij')
    #
    exe_tmp_wf = np.zeros((nstates, nspinorr, nspinorr, min(nk,block_size) ,
                           fft_box[0], fft_box[1], fft_box[2]),dtype=np.complex64)
    #
    exp_tmp_kL = np.zeros((min(nk,block_size) ,supercell[0],
                           supercell[1],supercell[2]),dtype=np.complex64)

    nblks = nk//block_size 
    nrem = nk%block_size
    if nrem > 0: nblks = nblks+1
    #
    pbar = tqdm(total=nk, desc="Ex-wf")
    #
    for ibk in range(nblks):
        ikstart = ibk*block_size
        ikstop = min(ikstart + block_size,nk)
        for ik in range(ikstart,ikstop):
            ## First get the electronic wfcs
            ik_ibz = kpt_idx[ik]
            isym = sym_idx[ik]
            wfc_tmp, gvecs_tmp = wfcdb.get_iBZ_wf(ik_ibz)
            wfc_tmp = wfc_tmp[:,elec_bnds[0]:elec_bnds[1],:,gvecs_iBZ_idx[ik_ibz]]
            gvecs_tmp = gvecs_tmp[gvecs_iBZ_idx[ik_ibz]]
            kvec = wfcdb.get_iBZ_kpt(ik_ibz)
            # get the rotated wf
            if isym != 0:
                sym_mat = wfcdb.ydb.sym_car[isym]
                time_rev = (isym >= len(wfcdb.ydb.sym_car
                                        ) / (1 + int(np.rint(wfcdb.ydb.time_rev))))
                kvec, wfc_tmp, gvecs_tmp = wfcdb.apply_symm(
                    kvec, wfc_tmp, gvecs_tmp, time_rev, sym_mat)

            kelec = kvec
            wfc_elec = wfc_tmp
            gvecs_elec = gvecs_tmp

            ## Do the same and get hole wfc
            ikhole = find_kpt(ktree, kelec-Qpt)
            ik_ibz = kpt_idx[ikhole]
            isym = sym_idx[ikhole]
            wfc_tmp, gvecs_tmp = wfcdb.get_iBZ_wf(ik_ibz)
            wfc_tmp = wfc_tmp[:,hole_bnds[0]:hole_bnds[1],:,gvecs_iBZ_idx[ik_ibz]]
            gvecs_tmp = gvecs_tmp[gvecs_iBZ_idx[ik_ibz]]
            kvec = wfcdb.get_iBZ_kpt(ik_ibz)
            # get the rotated wf
            if isym != 0:
                sym_mat = wfcdb.ydb.sym_car[isym]
                time_rev = (isym >= len(wfcdb.ydb.sym_car) / (1 + int(np.rint(wfcdb.ydb.time_rev))))
                kvec, wfc_tmp, gvecs_tmp = wfcdb.apply_symm(kvec, wfc_tmp, gvecs_tmp, time_rev, sym_mat)
            #
            khole = -kvec
            wfc_hole = wfc_tmp.conj()
            gvecs_hole = -gvecs_tmp

            if fix_particle == 'h':
                fx_kvec = khole
                fx_wfc = wfc_hole
                fx_gvec = gvecs_hole
                #
                ft_kvec = kelec
                ft_wfc = wfc_elec
                ft_gvec = gvecs_elec
            else :
                ft_kvec = khole
                ft_wfc = wfc_hole
                ft_gvec = gvecs_hole
                #
                fx_kvec = kelec
                fx_wfc = wfc_elec
                fx_gvec = gvecs_elec
            # compute
            ## NM : Donot perform FFT as we only need it for one point.
            exp_fx = np.exp(2*np.pi*1j*((fx_gvec + fx_kvec[None,:])@fixed_postion))
            fx_wfc *= exp_fx[None,None,None,:]
            fx_wfc = np.sum(fx_wfc,axis=-1) #(spin,bnd,spinor)
            ## for now only nspin = 1 works.
            ns1, nbndc, nspinorr, ng = ft_wfc.shape
            #if ft_ikpt not in prev_ikpts:
            ft_wfcr = wfcdb.to_real_space(ft_wfc.reshape(-1,nspinorr,ng),ft_gvec, grid=fft_box)
            ft_wfcr = ft_wfcr.reshape(ns1,nbndc,nspinorr,fft_box[0],fft_box[1],fft_box[2])
            exp_kx_r = np.exp(2*np.pi*1j*FFFboxs.reshape(-1,3)@ft_kvec)
            #
            if fix_particle == 'h':
                np.einsum('ncv,vy,cxijk->nxyijk',Akcv[:,ik,...],fx_wfc[0],ft_wfcr[0],
                          optimize=True,out=exe_tmp_wf[:,:,:,ik-ikstart])
            else :
                np.einsum('ncv,cx,vyijk->nxyijk',Akcv[:,ik,...],fx_wfc[0],ft_wfcr[0],
                          optimize=True,out=exe_tmp_wf[:,:,:,ik-ikstart])
            #
            exe_tmp_wf[:,:,:,ik-ikstart] *= exp_kx_r[...].reshape(FFFboxs.shape[:3])[None,None,None]
            exp_tmp_kL[ik-ikstart] = np.exp(1j*2*np.pi*np.einsum('...x,x->...',Lsupercells,ft_kvec))
            #
            # update progess bar
            pbar.update(1)
            #
        ## perform gemm operation 
        total_gemms_t = nstates*nspinorr**2
        exp_tmp_kL_tmp = exp_tmp_kL.reshape(len(exp_tmp_kL),-1)[:(ikstop-ikstart)].T
        exe_tmp_wf_tmp = exe_tmp_wf.reshape(nstates,nspinorr,nspinorr,-1,np.prod(fft_box))
        exe_tmp_wf_tmp = exe_tmp_wf_tmp[...,:(ikstop-ikstart),:]
        #
        for igemms in range(total_gemms_t):
            ii, jj, kk = np.unravel_index(igemms, (nstates,nspinorr,nspinorr))
            Ctmp = (exp_tmp_kL_tmp @ exe_tmp_wf_tmp[ii, jj, kk ])
            Ctmp = Ctmp.reshape(supercell[0], supercell[1], supercell[2],
                                fft_box[0], fft_box[1], fft_box[2])
            exe_wfc_real[ii, jj, kk ] += Ctmp.transpose(0,3,1,4,2,5)
    #
    exe_wfc_real   = exe_wfc_real.reshape(nstates,nspinorr,nspinorr,
                                            supercell[0]*fft_box[0],
                                            supercell[1]*fft_box[1],
                                            supercell[2]*fft_box[2])
    exe_wfc_real *= (1.0/np.prod(supercell))
    #
    # compute postioon of fixed particle in cart units 
    fixed_postion_cc = lat_vec@fixed_postion
    Lsupercells = Lsupercells.reshape(-1,3)#/np.array(supercell)[None,:]
    Lsupercells = Lsupercells@lat_vec.T
    atom_pos = Lsupercells[:,None,:] + wfcdb.ydb.car_atomic_positions[None,:,:]
    atom_pos = np.append(atom_pos.reshape(-1,3),fixed_postion_cc[None,:],axis=0)
    supercell_latvecs = lat_vec*np.array(supercell)[None,:]
    ## Make atomic numbers 
    atom_nums = np.zeros(len(atom_pos),dtype=wfcdb.ydb.atomic_numbers.dtype)
    attmp = atom_nums[:-1].reshape(-1,len(wfcdb.ydb.atomic_numbers))
    attmp[...]= wfcdb.ydb.atomic_numbers[None,:]
    atom_nums[-1] = 200
    #
    return supercell_latvecs,atom_nums,atom_pos,exe_wfc_real

#def compute_exc_wfc_real(path='.', bse_dir='SAVE', iqpt=1, nstates=[1],
#                          fixed_postion=[0,0,0], fix_particle='h', aveg=True, supercell=[1,1,1],
#                         wfcCutoffRy=-1, phase=False, block_size=256):
#    #
#    lattice = YamboLatticeDB.from_db_file(os.path.join(path, 'SAVE', 'ns.db1'))
#    filename = 'ndb.BS_diago_Q%d' % (iqpt)
#    excdb = YamboExcitonDB.from_db_file(lattice, filename=filename,
#                                                 folder=os.path.join(path, bse_dir),
#                                                 Load_WF=True, neigs=max(nstates))
#    # Load the wavefunction database
#    wfdb = YamboWFDB(path=path, latdb=lattice,
#                      bands_range=[np.min(excdb.table[:, 1]) - 1,
#                        np.max(excdb.table[:, 2])])
#    #
#    Akcv = excdb.get_Akcv()[min(nstates)-1:max(nstates)]
#    excQpt = excdb.car_qpoint
#    #
#    # Convert the q-point to crystal coordinates
#    Qpt = wfdb.ydb.lat @ excQpt

#    sc_latvecs, atom_nums, atom_pos, real_wfc = ex_wf2Real(Akcv, Qpt, wfdb, [np.min(excdb.table[:, 1]),
#                np.max(excdb.table[:, 2])], fixed_postion=fixed_postion,
#                          fix_particle=fix_particle, supercell=supercell, 
#                          wfcCutoffRy=wfcCutoffRy, block_size=block_size)
#    #
#    #
#    nstates_range = np.arange(nstates[0],nstates[1],dtype=int)
#    density = np.abs(real_wfc)**2

#    if fix_particle == 'h': name_file = 'electron'
#    else: name_file = 'hole'
#    if real_wfc.shape[1] != 1:
#        print("phase plot only works for nspin = 1 and nspinor == 1")
#        phase = False
#    if phase: 
#        phase = np.sign(real_wfc.real) #np.sign(np.angle(real_wfc))
#        density *= phase
#    if aveg:
#        real_wfc = np.sum(density,axis=(0,1,2))
#        real_wfc = real_wfc/np.max(np.abs(real_wfc))
#        write_cube('exe_wf_avg_%s_%d-%d.cube' %(name_file,nstates[0],nstates[1]),
#                   real_wfc, sc_latvecs, atom_pos, atom_nums, header='Real space exciton wavefunction')
#    else:
#        real_wfc = np.sum(density,axis=(1,2))
#        for i in range(len(real_wfc)):
#            real_wfc1 = real_wfc[i]/np.max(np.abs(real_wfc[i]))
#            write_cube('exe_wf_%s_%d.cube' %(name_file,nstates[i]), real_wfc1, sc_latvecs,
#                       atom_pos, atom_nums, header='Real space exciton wavefunction')



##if __name__ == "__main__":

