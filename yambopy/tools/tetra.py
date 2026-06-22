#
# Authors: FP
#
"""
Tetrahedron interpolation method for the k-integrals in DOS-like and 
spectral function-like expressions.

1. Based on
2. Inputs and outputs for tetrahedra construction
3. Inputs and outputs for integral calculation
4. Linear interpolation options
5. Analytical expressions for the k-integrals
6. Implementation
7. Usage example

:: Based on: 

    - "Optimized" tetrahedron method by [Kawamura et al. 
      Phys. Rev. B 89, 094515 (2014)] and its implementation 
      in Quantum ESPRESSO (`PW/src/tetra.f90`)
    - Integration formulas from
      [MacDonald et al. J. Phys. C: Solid State Phys. 12 (1979)]

:: Input and outputs for tetrahedra construction: 
                    
    - Inputs: check `get_tetrahedra_mesh` function for info
    - Outputs: you obtain a OptimizedTetrahedra object

:: Input and outputs for integral calculation: 

    - Inputs: check `spectra_tetrahedron` function for info
    - Outputs: you obtain the values of your integral expression

:: Linear interpolation options

    E_nk      : energies
    F_nk      : matrix elements
    i=1,2,3,4 : tetrahedra vertices

Linear (simple) case:

    E_n(k) -> E_n(\sum_i a_i k_i) = \sum_i a_i E_n(k_i)
    F_n(k) -> F_n(\sum_i a_i k_i) = \sum_i a_i F_n(k_i)

Optimized method (PUB):

    Here E_n(k_i) and F_n(k_i) are not exactly the values
    at the vertices of the tetrahedra, but the results of a
    polynomial interpolation linearized via least-square fits,
    producing a series of weights depending on polynomial
    degree. Everything is recast into a linear 
    transformation for E_nk and F_nk via the weights matrix.

:: Analytical expressions for the k-integrals

    The expressions for D_T and A_T can be found in PAPER, EQUATIONS

    T         : tetrahedron index
    i=1,2,3,4 : vertices of tetrahedron T

DOS case:

    D(w) = 1/Nk \sum_n \int_dk \delta(w-E_nk)

    becomes

    D(w) = 1/NT \sum_n \sum_T D_T(w,w1,w2,w3,w4)

Spectral function case:

    A(w) = 1/Nk \sum_n \int_dk F_nk \delta(w-E_nk)

    becomes

    A(w) = 1/NT \sum_n \sum_T D_T(w,w1,w2,w3,w4) \sum_i A_T^i(w,w1,w2,w3,w4) F^i_T

:: Implementation

Conventions:

    We construct the tetrahedra using the same ordering 
    convention of the QE implementation by M. Kawamura (PATH).

    In this way, the DOS spectra ought to exactly correspond 
    to QE results at any value of Nk.

    Other conventions are possible, in which case DOS spectra 
    ought to correspond to QE results only in the converged 
    limit of large Nk.

Parallelization:
    
    The spectra evaluation is parallelized with `joblib`, 
    by default over w values, otherwise over n loop.

:: Usage example
```
import numpy as np
from yambopy.tools.tetra import *
from yambopy import YamboLatticeDB   # (optional,recommended) to get `red_kpoints` and `rlat`
from yambopy import YamboElectronsDB # (optional,use-case) Energies for DOS
from qepy import ProjwfcXML          # (optional,use-case) Projections for PDOS
from yambopy import YamboWFDB        # (optional,use-case) Projections for PDOS (BZ expansion)

# Parallel jobs
njobs=8

# Database info
ns1 = 'SAVE/ns.db1'
prefix = 'hBN'
pdos_outfile = 'pdos.out'

#
# DOS case
#
# You need: (i) reciprocal space geometry, (ii) poles of DOS
#

nk1,nk2,nk3 = [6,6,2]
lat = YamboLatticeDB.from_db_file(filename=ns1) # (i) reciprocal space geometry
el  = YamboElectronsDB.from_db_file(filename=ns1)
energies=el.eigenvalues_ibz[0][lat.kpoints_indexes] # (ii) Poles of SF / DOS

# Initialize energies
emin=2.
emax=14.
estep=0.02
w = np.arange(emin,emax+estep,estep)
nfreqs = len(w)

# Run tetrahedron method (optimized tetra)
## Build tetrahedra connections (fast)
tetra = get_tetrahedra_mesh(nk1,nk2,nk3,lat.red_kpoints,lat.rlat*lat.alat[0])
## Calculate integral (parallel calc)
DOS   = spectra_tetrahedron(energies,tetra,w,njobs=njobs)

#
# Spectral Function case (PDOS in this specific example)
#
# You need: (i) rec. space geometry, (ii) poles of SF, (iii) matrix elements of SF
#

# Preparation of expanded projections
projwfc = ProjwfcXML(prefix,output_filename=pdos_outfile)
wfdb    = YamboWFDB(path='.', latdb=lat)
proj_BZ = np.abs( projwfc.rotate_proj(wfdb) )**2.

# Selection of PDOS
p_states = []
for ip in range(projwfc.nproj):
    if projwfc.states[ip]['l']==1: p_states.append(ip)
p_proj_BZ = np.sum(proj_BZ[:,p_states,:],axis=1)   # (iii) Matrix elements of SF

# Run tetrahedron method (optimized + interpolation of matrix el. [default])
p_PDOS=spectra_tetrahedron(energies,tetra,w,matels=p_proj_BZ,njobs=njobs)
# Run tetrahedron method (linear + average of matrix el.)
#p_PDOS=spectra_tetrahedron(energies,tetra,w,matels=p_proj_BZ,linear=lin,interp_matels=False,njobs=njobs)
```

# You can plot (w,DOS) and (w,p_PDOS)
"""
import numpy as np
from yambopy.tools.degeneracy_finder import find_degeneracy_evs
from yambopy.tools.citations import citation
from yambopy.kpoints import build_ktree,find_kpt
from joblib import Parallel,delayed
from tqdm import tqdm

class OptimizedTetrahedra:
    """
    This class contains the quantities needed to apply
    the optimized tetrahedron method for interpolation of
    reciprocal-space integrals

    :: tetra -> np array (ntetra,20) of k-point indices.
                For each tetrahedron, it contains the 20
                corresponding indices (first 4: the actual vertices
                of the tetrahedron, then: the remaining 
                additional points of the optimized method)

    :: w_sp  -> np array (4, 20). The values of the weights
                used in the optimized method to reconstruct
                effective energies from more than the 4 tetra
                vertices

    :: ntetra  -> number of tetrahedra: 6*nkpoints_bz
    :: nntetra -> number of total kpoints per tetrahedron (4+16).
                  Used in the optimized method, hard-set to 20.

    :: verts   -> [optional] np array (6, 20, 3). 
                  If supplied, the basic tetrahedra vertex structure
                  not assigned to k-point indices
    """
    def __init__(self,tetra,w_sp,ntetra,nntetra=20,verts=None):
        self.tetra   = tetra
        self.w_sp    = w_sp
        self.ntetra  = ntetra
        self.nntetra = nntetra
        if verts is not None: self.verts = verts

@citation("M. Kawamura et al. Phys. Rev. B 89, 094515 (2014)")
def get_tetrahedra_mesh(nk1,nk2,nk3,red_kpts,rlat_cc):
    """
    This function returns an OptimizedTetrahedra object for
    DOS-like and spectral-function-like interpolations with 
    the optimized tetrahedron method.

    Parameters
    ----------
    * nk1,nk2,nk3 -> Monkhorst-Pack grid size
    * red_kpts    -> kpts in FULL BZ in fractional coords.
    * rlat_cc     -> reciprocal lattice vectors in QE cart. coord.: 
                                lat.rlat*lat.alat[0]
                     (used only for QE-compliant tetra initialization)

    Output
    -------
    * OptimizedTetrahedra object
    
    It constructs the tetrahedra in which k-space is subdivided
    (6 per "cube" volume dk), getting the location of their 4 
    vertices in k-indices of the full BZ. 

    It gets also the additional 16 points (total 20 pts) that are
    used in the optimized method, as well as the corresponding
    special weights for effective energy and mat. el. evaluation.

    Check comments to understand what is going on.
    """

    # We start from the unit "cube" in reciprocal space,
    # connecting 8 closest k-points (the vertices) of
    # the Monkhorst-Pack grid in cart. coord.
    d3k = np.array([ rlat_cc[0]/nk1 , rlat_cc[1]/nk2, rlat_cc[2]/nk3 ])

    # We want to split the cube into tetrahedra sharing a diagonal. 
    # As a starting point, we identify the shortest diagonal
    d3k_diags = np.array([
        -d3k[0] + d3k[1] + d3k[2],
         d3k[0] - d3k[1] + d3k[2],
         d3k[0] + d3k[1] - d3k[2],
         d3k[0] + d3k[1] + d3k[2] ])

    len_diag   = np.sum(d3k_diags**2, axis=1)
    i_min_diag = np.argmin(len_diag)

    # Now we build the tetrahedra by finding their 4 vertices
    verts0 = np.zeros(4,dtype=int)
    verts0[i_min_diag] = 1 # first vertex

    # 6 tetrahedra per cube
    # Optimized method: 20 additional points per tetrahedron (see later)
    verts = np.zeros((6,20,3), dtype=int)

    # Directional steps to find other vertices
    directions = np.eye(4,dtype=int)
    directions[i_min_diag,i_min_diag] = -1 # first direction

    # If we start from a cube vertex v0 and we move along cube
    # edges ei, we find 3 other vertices: this is one tetrahedron
    # v0          -> first vertex
    # v0+e1       -> second vertex
    # v0+e1+e2    -> third vertex
    # v0+e1+e2+e3 -> fourth vertex
    # So we get the vertices of the 6 tetrahedra in the cube
    i_verts = 0
    for ix in range(3):
        for iy in range(3):
            if ix==iy: continue
            for iz in range(3):
                if iz==ix or iz==iy: continue

                verts[i_verts,0]=verts0[:3] # e.g. [1,0,0]: vertex
                verts[i_verts,1]=verts[i_verts,0]+directions[:3,ix]
                verts[i_verts,2]=verts[i_verts,1]+directions[:3,iy]
                verts[i_verts,3]=verts[i_verts,2]+directions[:3,iz]

                i_verts+=1

    # Now there are other 16 additional points/vertices that 
    # correspond to the optimized method by Kawamura et al
    # For example, points external to the tetrahedron are stored,
    # such as v4 = v0 + (v0-v1), and many others
    verts[:, 4]  = 2 * verts[:, 0] - verts[:, 1]
    verts[:, 5]  = 2 * verts[:, 1] - verts[:, 2]
    verts[:, 6]  = 2 * verts[:, 2] - verts[:, 3]
    verts[:, 7]  = 2 * verts[:, 3] - verts[:, 0]
    verts[:, 8]  = 2 * verts[:, 0] - verts[:, 2]
    verts[:, 9]  = 2 * verts[:, 1] - verts[:, 3]
    verts[:, 10] = 2 * verts[:, 2] - verts[:, 0]
    verts[:, 11] = 2 * verts[:, 3] - verts[:, 1]
    verts[:, 12] = 2 * verts[:, 0] - verts[:, 3]
    verts[:, 13] = 2 * verts[:, 1] - verts[:, 0]
    verts[:, 14] = 2 * verts[:, 2] - verts[:, 1]
    verts[:, 15] = 2 * verts[:, 3] - verts[:, 2]
    verts[:, 16] =     verts[:, 3] - verts[:, 0] + verts[:, 1]
    verts[:, 17] =     verts[:, 0] - verts[:, 1] + verts[:, 2]
    verts[:, 18] =     verts[:, 1] - verts[:, 2] + verts[:, 3]
    verts[:, 19] =     verts[:, 2] - verts[:, 3] + verts[:, 0]    

    # So now the idea is that all these new points are used to
    # construct a linear fit of polynomially interpolated energies, 
    # to get new energies with values 
    # that do not match exactly the vertices of the six original 
    # tetrahedra in k-space: this is to correct the errors 
    # in the standard linear interpolation. 
    # In order to find these effective
    # energies, we need to construct special weights w_sp.
    #
    # If i are the 4 tetra vertices, and k runs over the 20 special
    # points, we have: 
    #                  E_i = sum_k w_ik E0_k
    # Where E0_k is the actual calculated energy at that k-point, 
    # and E_i is the effective energy.
    w_sp = np.zeros((4,20), dtype=float)

    # Values of the weights are taken from Kawamura's implementation:
    # They are fixed by the choice of a third-order polynomial for
    # the coefficient expansion (a linear least-square fit is then
    # performed, leading to the Eq. above).
    w_sp[:, 0:4] = np.array([ [1440,    0,   30,    0],
                              [   0, 1440,    0,   30],
                              [  30,    0, 1440,    0],
                              [   0,   30,    0, 1440] ])

    w_sp[:, 4:8] = np.array([ [-38,   7,  17, -28],
                              [-28, -38,   7,  17],
                              [ 17, -28, -38,   7],
                              [  7,  17, -28, -38] ])

    w_sp[:, 8:12] = np.array([ [-56,   9, -46,   9], 
                               [  9, -56,   9, -46], 
                               [-46,   9, -56,   9],
                               [  9, -46,   9, -56] ])

    w_sp[:, 12:16] = np.array([ [-38, -28,  17,   7],
                                [  7, -38, -28,  17],
                                [ 17,   7, -38, -28],
                                [-28,  17,   7, -38] ])

    w_sp[:, 16:20] = np.array([ [-18, -18,  12, -18],
                                [-18, -18, -18,  12],
                                [ 12, -18, -18, -18],
                                [-18,  12, -18, -18] ])

    w_sp /= 1260.0

    # Now that we have everything set, we need to map the tetra
    # vertices to the actual k-point indices, 
    # in order to get their energies!
    nkpoints = nk1*nk2*nk3
    ntetra   = 6*nkpoints
    Nk       = np.array([nk1,nk2,nk3])

    # For each i_tetra, this will contain all the 20 real k-indices
    tetra = np.zeros((ntetra,20), dtype=int)

    # Matching tree
    ktree = build_ktree(red_kpts)

    # Assign indices to all tetrahedra 4+16 vertices
    i_all_tetra = 0
    for kpoint in red_kpts: # for each point
        for i_tetra in range(6): # for each tetrahedron

            # Add the 20 points as directions in fractional coordinates
            kpts_T = kpoint + verts[i_tetra]/Nk
            # Map to the BZ indices
            ik_BZ_tetra = [ find_kpt(ktree,kpts_T[i])  for i in range(20) ]
            # Finally, the 20 full BZ indices associated with i_tetra
            tetra[i_all_tetra] = ik_BZ_tetra

            i_all_tetra += 1

    # This small object contains everything necessary
    return OptimizedTetrahedra(tetra=tetra,w_sp=w_sp,ntetra=ntetra,verts=verts)

@citation("M. Kawamura et al. Phys. Rev. B 89, 094515 (2014)")
@citation("A. H. MacDonald et al. J. Phys. C: Solid State Phys. 12 (1979)")
def spectra_tetrahedron(energies,tetra_mesh,freqs,nspin=2,matels=None,interp_matels=True,linear=False,njobs=1,par_mode='auto'):
    """
    Calculate a k-integral with the tetrahedron method:

    D(w;q) = \sum_n \int_dk |C_nk(q)|^2 \delta( w-E_nk(q) )

    If the expression depends on multiple state and momenta 
    indices and energies (e.g., lifetimes, spectral functions 
    at finite momentum), you can provide a generalized mapping 
    to bring back the expression in the above shape with
    `n` generalized state index and `E_nk` generalised pole energy.

    This function prepares the caculation and sets the 
    parallelization with `joblib`, then calls 
    `calculate_spectra_tetrahedron` for each w value in `freqs`
    
    Parameters
    -----------
    energies   -> pole energies in eV [e_nk], 
                  must be an array of shape (N_momenta, N_states)
    freqs      -> Evaluation energies in eV [w] (np array)
    matels     -> matrix elements squared [|C_nk|^2, optional], 
                  must be an array of shape (N_momenta,N_states).
                  * If None:     DOS calculation.
                  * If not None: SPECTRAL FUNCTION calculation.
    interp_matels -> Only used for Spectral Function calculation.
                  * If True (default): linear interpolation of |C_nk|^2
                  * If False: take average value of |C_nk|^2 
                              inside each tetrahedron 
    tetra_mesh -> OptimizedTetrahedra object
    linear     -> Use simple tetrahedron method instead of 
                  optimized [default: False]
    njobs      -> no. of parallel jobs (requires `joblib`). 
                  [default: serial]
    par_mode   -> which parallel loop: 
                  'auto' (default), 
                  'w' (energies), 
                  'b' (states).
                  * If njobs > 1 and nfreqs>=nstates : 'auto'=='w'
                  * If njobs > 1 and nfreqs<nstates  : 'auto'=='b'
    
    Outputs
    -----------
    Values of D( w ) as np array of len(freqs)

    """

    # Check if energies are array or just one value
    try: nfreqs = len(freqs)
    except TypeError:
        freqs = np.array([freqs])
        nfreqs = len(freqs)
     
    # System size
    nstates  = energies.shape[1]
    nmomenta = energies.shape[0]
    ntetra   = tetra_mesh.ntetra

    # Parallel options
    if par_mode not in ['auto','w','b']: par_mode='auto'
    if par_mode=='auto' and nfreqs>=nstates: par_mode='w'
    if par_mode=='auto' and nfreqs<nstates: par_mode='b'
    # Parallel checks
    if par_mode=='w' and njobs>nfreqs:
        return ValueError(f'[ERROR] {njobs} jobs requested but only {nfreqs} energy steps')
    if par_mode=='b' and njobs>nstates:
        return ValueError(f'[ERROR] {njobs} jobs requested but only {nstates} eigenvalues')

    # Print info
    if linear: tmode='linear'
    else:      tmode='improved'
    if matels is None: ftype = 'DOS'
    else:              ftype = 'Spectral Function'
    print(f":: Tetrahedron method ({tmode}) ::")
    if not interp_matels: print(f"   {ftype} calculation (averaged matrix elements)")
    else:                 print(f"   {ftype} calculation")
    print(f"   Energy steps: {nfreqs}, states: {nstates}, momenta: {nmomenta}, tetrahedra: {ntetra}")
    if njobs>1 and par_mode=='w':
        print(f"   Parallelization over energy steps ({njobs} jobs)")
    if njobs>1 and par_mode=='b':
        print(f"   Parallelization over states ({njobs} jobs)")

    # Run calculation
    if par_mode=='w':
        # Parallelization of frequency steps (suggested in spectral function case)
        S =  np.array( list( tqdm( Parallel(return_as="generator",n_jobs=njobs,backend='loky')(delayed(calculate_spectra_tetra)(energies,tetra_mesh,freqs[iw],nspin=nspin,matels=matels,interp_matels=interp_matels,linear=linear,njobs=1) for iw in range(nfreqs)), total=nfreqs, desc="Tetrahedron interpolation" ) ) )

    if par_mode=='b':
        S = np.zeros(nfreqs)
        # Parallelization over eigenvalues (suggested in lifetimes case with only one freq)
        for iw in range(nfreqs):
            S[iw] = calculate_spectra_tetra(energies,tetra_mesh,freqs[iw],nspin=nspin,matels=matels,interp_matels=interp_matels,linear=linear,njobs=njobs,iw=iw)

    return S

def calculate_spectra_tetra(energies,tetra_mesh,w,nspin=2,matels=None,interp_matels=True,linear=False,njobs=1,iw=None):
    """
    Calculate DOS/spectral function/lifetime using the tetrahedron method at
    the evaluation energy `w`.

    Uses tetrahedron formulas. Both linear and optimised methods.

    Actual calculation per state in internal function `k_integral_at_k_and_n`.
    
    Parameters
    -----------
    energies   -> pole energies in eV, must be an array of 
                  shape (N_momenta, N_states)
    w          -> Evaluation energy in eV
    nspin      -> Spin factor depending on system 
    matels     -> matrix elements squared, must be an array
                  of shape (N_momenta,N_states)
    interp_matels -> whether to use full spectral function expression or <M>_T*DOS_T
    tetra_mesh -> OptimizedTetrahedra object
    linear     -> Use simple tetrahedron method instead of optimized [default: False]
    njobs      -> Parallelize the states loop [default: 1 (serial)]
    iw         -> energy step index (used only if parallel calc with nfreqs>1)

    Output
    -----------
    Value of spectrum at w (float scalar)
    
    """
    def k_integral_at_k_and_n(istate):
        """
        Function that performs the actual calculation, compatible with `joblib` call.

        - Obtain F_n(w=w0) = int_dk G_nk(w0) : `istate` is `n` index
        - Full result will be F(w) = sum_n F_n(w) for all w
        """
        if linear: # only the 4 tetra vertices
            E_eff = energies[tetra[:,:4],istate]
            if spectral_function: M_eff = matels[tetra[:,:4],istate]
        else: # the 20 vertices + effective energies
            # The 20 energy values for each optimized tetrahedron per band
            E0 = energies[tetra,istate]
            # Now we calculate the actual effective energies as
            # E_i = sum_ik w_ik E0_k
            E_eff = np.einsum('ik,tk->ti',opt_weights,E0)
            #Same thing for the matrix elements!
            if spectral_function:
                M0 = matels[tetra,istate]
                M_eff = np.einsum('ik,tk->ti',opt_weights,M0)

        # Sort vertex energies per tetrahedron, e.g.:
        # E1<E2<E3<E4
        E_order = np.argsort(E_eff, axis=1)
        # Reorder energies...
        E_eff = np.take_along_axis(E_eff, E_order, axis=1)
        # ... and of course reorder mat elements
        if spectral_function: M_eff = np.take_along_axis(M_eff, E_order, axis=1)
        
        # Now we have the effective energy at each vertex per tetra:
        E1, E2, E3, E4 = E_eff.T
       
        # And the effective matrix elements:
        if spectral_function:   
            if interp_matels: M1, M2, M3, M4 = M_eff.T # linear interpolation
            else:             M_mean = np.mean(M_eff, axis=1) # average per tetra

        # In the tetrahedron method, the momentum integral is done
        # analytically. The integral evaluation is split 
        # into the regions R1: E1<w<E2, R2: E2<=w<E3, R3: E3<=w<E4.

        # For each tetrahedron, indices belonging to the various regions
        R1 = (w >  E1) & (w < E2)
        R2 = (w >= E2) & (w < E3)
        R3 = (w >= E3) & (w < E4)
        #R4 = (w >= E4 ) # zero for DOS and spectral functions
        #R5 = (w <= E1 ) # zero for DOS and spectral functions

        # Preparation for analytical formulas
        with np.errstate(divide='ignore', invalid='ignore'):
            f21 = ( w - E1 ) / ( E2 - E1 )
            f31 = ( w - E1 ) / ( E3 - E1 )
            f41 = ( w - E1 ) / ( E4 - E1 )
            f12 = 1. - f21
            f32 = ( w - E2 ) / ( E3 - E2 )
            f42 = ( w - E2 ) / ( E4 - E2 )
            f13 = 1. - f31 
            f23 = 1. - f32
            f43 = ( w - E3 ) / ( E4 - E3 )                         
            f14 = 1. - f41
            f24 = 1. - f42
            f34 = 1. - f43

        # Now the analytical formulas
        DOS = np.zeros(ntetra)

        # DOS contributions
        DOS[R1] = 3.* f21[R1]*f31[R1]*f41[R1] / (w-E1[R1]) 
        DOS[R2] = 3.* (f23[R2]*f31[R2] + f32[R2]*f24[R2]) / (E4[R2]-E1[R2])
        DOS[R3] = 3.* f14[R3]*f24[R3]*f34[R3] / (E4[R3]-w)

        # DOS
        if not spectral_function:
            # This is the result of the k-integration, now summed over
            # k-regions (np.sum): \int_dk \delta(w-Ek)
            return np.sum( DOS )

        # Spectral function (matrix element contribution)
        elif not interp_matels:
            # M_mean taken outside of the integral for each tetra, then:
            return np.sum( M_mean * DOS )

        else:
            # Full formula with linear interpolation of M inside tetra
            V1 = np.zeros(ntetra) # vertex 1
            V2 = np.zeros(ntetra) # vertex 2
            V3 = np.zeros(ntetra) # vertex 3
            V4 = np.zeros(ntetra) # vertex 4
            
            # R1
            V1[R1] = ( f12[R1]+f13[R1]+f14[R1] )/3.
            V2[R1] = f21[R1]/3.
            V3[R1] = f31[R1]/3.
            V4[R1] = f41[R1]/3.

            # R2
            DE = DOS[R2] * ( E4[R2]-E1[R2] ) # cumbersome denominator
            V1[R2] = f14[R2]/3. + ( f13[R2]*f31[R2]*f23[R2] ) / DE
            V2[R2] = f23[R2]/3. + ( f24[R2]**2. *f32[R2] )    / DE
            V3[R2] = f32[R2]/3. + ( f31[R2]**2. *f23[R2])     / DE
            V4[R2] = f41[R2]/3. + ( f42[R2]*f24[R2]*f32[R2] ) / DE

            # R3
            V1[R3] = f14[R3]/3.
            V2[R3] = f24[R3]/3.
            V3[R3] = f34[R3]/3.
            V4[R3] = ( f41[R3]+f42[R3]+f43[R3] )/3.
          
            # This is the expression of the spectral function
            # after the integral is performed including the linear
            # variations of the matrix elements: \int_dk Mk * \delta(w-Ek)
            # Sum over all k-regions R1,R2,R3 -> np.sum
            return np.sum( (V1*M1 + V2*M2 + V3*M3 + V4*M4) * DOS )

    # Preparation
    tetra       = tetra_mesh.tetra
    opt_weights = tetra_mesh.w_sp
    ntetra      = tetra_mesh.ntetra

    nmomenta = energies.shape[0]
    nstates  = energies.shape[1]

    # Tetrahedra in full BZ, values must match
    assert nmomenta == ntetra/6

    if matels is None: 
        spectral_function = False
    else:
        spectral_function = True
        assert np.isrealobj(matels)
        assert matels.shape[0] == ntetra/6
        assert matels.shape[1] == nstates
        # As pre-interpolation step, we average over mat. elements
        # corresponding to degenerate states
        # FP: so far this cannot be toggled by the user, but it may have negligible effects
        matels = deg_average(matels,energies)    
    
    # Value of spectrum at energy w
    spectrum_w = 0.0
    
    if njobs==1: # Avoid joblib function if we may be using it for energy steps
        for istate in range(nstates): 
            # Sum over band states (+=)
            spectrum_w += k_integral_at_k_and_n(istate)

    if njobs>1: # We are not parallelizing over energies so we go over states
            spectrum_w_n =  np.array( list( tqdm( Parallel(return_as="generator",n_jobs=njobs,backend='loky')(delayed(k_integral_at_k_and_n)(istate) for istate in range(nstates)), total=nstates, desc=f"Tetrahedron interpolation @freq{iw}" ) ) )
            # Sum over band states
            spectrum_w = np.sum(spectrum_w_n)

    return spectrum_w*nspin/ntetra

def deg_average(values,energies):
    """
    Average matrix elements over degenerate subspaces.

    * Step 1: find indices of degenerate states in `energies`
    * Step 2: average `values`
    """

    new_values = np.copy(values)
    for ik in range(len(energies)):
        # Step 1
        deg_at_k = find_degeneracy_evs(energies[ik]) # List of np arrays
        # Step 2
        if deg_at_k is not None:
            for subspace in deg_at_k:
                if len(subspace)<=1: continue
                v_avg = np.mean(values[ik,subspace])
                new_values[ik,subspace]=v_avg

    return new_values
