# Yambo Input variables


All the default input flags parsed to ``yambo`` 

## RunLevels
- chi: [R][CHI] Dyson equation for Chi.
- bse: [R][BSE]Bethe Salpeter Equation
- tddft: [R][K] Use TDDFT kernel
- cohsex:  [R][Xp] COlumb Hole Screened EXchange
- Xx:       [R][Xx] Dynamical Response Function
- em1s:    [R][Xs] Statically Screened Interaction
- em1d:    [R][X] Dynamically Screened Interaction
- ppa:     [R][Xp] Plasmon Pole Approximation for the Screened Interaction
- mpa:     [R][Xm] Multi Pole Approximation for the Screened Interaction
- el_el_corr    [R] Electron-Electron Correlation

## RT SCATT
- el_el_scatt: [R] Electron-Electron Scattering
- el_photon_scatt: [R] Electron-Photon   Scattering
- 
## RT SCATT + ELPH

- el_ph_scatt: [R] Electron-Phonon   Scattering

## RT SCATT + PHEL
- ph_el_scatt [R] Phonon-Electron Scattering

## ELPH

- el_ph_corr   [R] Electron-Phonon Correlation  
- BSEscatt:    [KPT] Compute extended k/q scatering
- ElPhRndNq: [ELPH] Read random Q-pointselph_nQ_used
- EkpqShFact. [ELPH] E(k+q) Interpolation shell factor (used only with double-grid)

## PHEL 
- ph_el_corr   [R] Phonon-Electron Correlation   


## QED 

- el_photon_corr: [R] Electron-Photon Correlation 
- QEDRLvcs: [QED] Vector-Potential G-vectors

### CPUs
#### GPL variables

- StdoHash: [IO] Live-timing Hashes
- MaxGvecs: [INI] Max number of G-vectors planned to useng_closed
- Gthresh: [INI] Accuracy for energy shell separation [-1. automatic]
- FFTGvecs: [FFT] Plane-waveswf_ng 
- NonPDirs: [X/BSS] Non periodic chartesian directions (X,Y,Z,XY...)
- MolPos:  [X/BSS] Molecule coord in supercell, 0.5 is the middlemolecule_position
- K_grids: [KPT] Select the grids (X=response, S=sigma, C=collisions, B=bse) [default X S C B]
- IkSigLim: [KPT] QP K-points indices
- IkXLim:  [KPT] X grid last k-point index
- Nelectro: Electrons numbernel
- ElecTemp: Electronic TemperatureTel
- OccTresh: Occupation treshold (metallic bands)
- BoseTemp: Bosonic TemperatureBose_Temp
- BoseCut [BOSE] Finite T Bose function
- WFbuffIO: [IO] Wave-functions buffered I/Overb_level=V_io
- NoDiagSC: New setup for non-diagonal supercells

#### Parallel Setups
##### OPENMP
- K_Threads       [OPENMP/BSK] Number of threads for response function
- X_Threads       [OPENMP/X] Number of threads for response functions
- DIP_Threads     [OPENMP/X] Number of threads for dipoles
- SE_Threads      [OPENMP/GW] Number of threads for self-energy
- RT_Threads      [OPENMP/RT] Number of threads for real-time
- NL_Threads      [OPENMP/NL] Number of threads for nl-optics

###### MPI
- NLogCPUs     [PARALLEL] Live-timing CPU`s (0 for all)n_log_CPUs
## I/O
- DBsIOoff:    [IO] Space-separated list of DB with NO I/O. DB=(DIP,X,HF,COLLs,J,GF,CARRIERs,OBS,W,SC,BS,ALL)
- DBsFRAGpm:   [IO] Space-separated list of +DB to FRAG and -DB to NOT FRAG
## S xc : related to exchange correlation functionals
- LifeTrCG: [GW] [o/o] Lifetime transition reduction
- HARRLvcs: [HA] Hartree     RL components
- EXXRLvcs: [XX] Exchange    RL components
- UseNLCC:  [XC] If present, add NLCC contributions to the charge density
- VXCRLvcs: [XC] XCpotential RL components
- CORRLvcs: [GW] Correlation RL components

### SC
- GbndRnge: [GW] G[W] bands range
- UseEbands: [GW] Force COHSEX to use empty bands
- ALLGexx: [XX] Force the use use all RL vectors for the exchange part
- ALLGHAR: [HA] Force the use use all RL vectors for the Hartree potential


- GbndRnge: [GW] G[W] bands range
- UseEbands: [GW] Force COHSEX to use empty bands
- GDamping: [GW] G[W] damping
- #else
- GDamping: [GW] G[W] damping
- GDmRnge:  [GW] G_gw damping range
- dScStep:  [GW] Energy step to evaluate Z factors
- DysSolver: [GW] Dyson Equation solver ("n","s","g","q")
- GEnSteps: [GW] Green`s Function (GF) energy steps
- GEnRnge:  [GW] GF energy range
- GEnMode:  [GW] GF energy mode ("centered","absolute"). "Centered" around the bare energy
- GTermKind: [GW] GW terminator ("none","BG" Bruneval-Gonze,"BRS" Berger-Reining-Sottile)
- GTermEn:  [GW] GW terminator energy (only for kind="BG")
- NewtDchk:  [GW] Test dSc/dw convergence
- ExtendOut: [GW] Print all variables in the output file
- QPsymmtrz: [GW] Force symmetrization of states with the same energy

#### *_Xx Xs Xd Xp Xm_*

- Chimod:       [X] IP/Hartree/ALDA/LRC/PF/BSfxc
- ChiLinAlgMod: [X] inversion/lin_sys,cpu/gpuChi_linalg_mode
- XTermKind: [X] X terminator ("none","BG" Bruneval-Gonze)
- XTermEn:  [X] X terminator energy (only for kind="BG")
- DrClassic: [X] Use a classical model for the drude term
- mpERdb:  [Xm] Write to disk MPA poles and residues
- WriteXo: [X] Write on the-fly the IP response function

#### *_BSE/BSK_*
- BSEmod:  [BSE] resonant/retarded/coupling
- Lkind:   [BSE,X] bar(default)/fullBSE
- BSEBands: [BSK] Bands range
- BSENGBlk: [BSK] Screened interaction block size [if -1 uses all the G-vectors of W(q,G,Gp)]
- BSENGexx: [BSK] Exchange components
- BSENGfxc: [BSK] Fxc components 
- FxcCutScheme [TDDFT] ("none","lower_Gmax","lower_GmGp","full_grid")
- BSEEhEny: [BSK] Electron-hole energy range
- BSKCut  [BSK] Cutoff on the BSE Kernel
- BSKIOmode: [BSK] ("1D_linear"/"2D_standard" + "norestart")
- BSKmod  [BSE] IP/Hartree/HF/ALDA/SEX/BSfxcBSK_mode
- Gauge  [BSE/X] Gauge (length|velocity)
- NoCondSumRule [BSE/X] Do not impose the conductivity sum rule in velocity gauge
- MetDamp       [BSE] Define $//slash//w+=sqrt(//slash//w*(//slash//w+i//slash//eta))$
- BSSmod    [BSS] (h)aydock/(d)iagonalization/(s)lepc/(i)nversion/(t)ddft
- BSEprop   [BSS] Can be any among abs/jdos/kerr/magn/dich/photolum/esrt
- BSEdips   [BSS] Can be "trace/none" or "xy/xz/yz" to define off-diagonal rotation plane
- BSSInvMode: [BSS] Inversion solver modality `(f)ull/(p)erturbative`
- BSSInvPFratio: [BSS] Inversion solver. Ratio between the number of frequencies solved pert/full
- BLongDir  [BSS] [cc] Electric Field
- BEnRange  [BSS] Energy range
- BDmRange  [BSS] Damping range
- BSHayTrs  [BSS] Relative [o/o] Haydock threshold
- BSHayItrMAX: [BSS] MaX number of Haydock iterations
- BSHayItrIO [BSS] Iterations between IO for Haydock restart
- BSEPSInvTrs [BSS EPS] Inversion treshold
- BSPLInvTrs  [BSS PL] Inversion treshold
- BEnSteps  [BSS] Energy steps
- DrudeWBS  [BSE] Drude plasmon
- FxcLibxc: [BSK] force computing Fxc via libxc
- WehDiag [BSK] diagonal (G-space) the eh interaction
- WehCpl  [BSK] eh interaction included also in coupling
- WRbsWF  [BSS] Write to disk excitonic the WFs
### SLEPC and NOT NL

- BSSSlepcApproach: [SLEPC] Approach ("Krylov-Schur","Generalized-Davidson","Jacob-Davidson")
- BSSNEig       [SLEPC] Number of eigenvalues to compute
- BSSEnTarget   [SLEPC] Target energy to find eigenvalues
- BSSSlepcMaxIt [SLEPC] Maximum number of iterations
- BSSSlepcPrecondition: [SLEPC] Precondition technique (none|preonly+jacobi|bcgs+jacobi)
- BSSSlepcExtraction: [SLEPC] Extraction technique (ritz|harmonic)
- BSSSlepcNCV   [SLEPC] Dimension of the subspace
- BSSSlepcTol   [SLEPC] Tolerance for the iterative solver
- BSSSlepcMatrix: [SLEPC] Store slepc matrix, faster but more memory consuming

### SC
- ALLGexx [BSS] Force the use use all RL vectors for the exchange part
- BSHayTer: [BSS] Terminate Haydock continuos fraction
- Reflectivity [BSS] Compute reflectivity at normal incidence
- BSSPertWidth [BSS] Include QPs lifetime in a perturbative way
- BSSInvKdiag: [BSS] In the inversion solver keep the diagonal kernel in place

#### *_Fxc_*
- FxcGRLc  [TDDFT] XC-kernel RL size
- LRC_alpha: [TDDFT] LRC alpha factor
- PF_alpha: [TDDFT] PF alpha factor approximation: CUR/EMP/RBO/JGMF
- LRC_beta [TDDFT] LRC beta factor

### Optics: Q=0 directions average

- OptDipAverage [Xd] Average Xd along the non-zero Electric Field directions

### Optics: large Q momenta

- Qdirection [Xd] Transferred momentum direction (iku)
- QShiftOrder: [Xd] Pick-up the (QShiftOrder)th q+G

### BSE properties

- AnHall  , [BSE] Add the anomalous Hall effect to eps if using length gauge
- PL_weights: [PL] [cc] Weights of the carthesian components of the emitted radiation

#### BSE: Real-Time

- RTOccMode: [RT-BSE] (K)ernel/(R)esiduals. BSE components to be corrected with the TD occupations
- ForceEqTrans: [RT-BSE] Use only equilibrium transitions

### Double Grid(s)

- DbGdQsize [X,DbGd][o/o] Percentual of the total DbGd transitions to be used

### RIM

- Em1Anys [RIM] X Y Z Static Inverse dielectric matrix Anysotropy
- IDEm1Ref: [RIM] Dielectric matrix reference component 1(x)/2(y)/3(z)
- RandQpts: [RIM] Number of random q-points in the BZ
- RandGvec: [RIM] Coulomb interaction RS componentsRIM
- QpgFull [F RIM] Coulomb interaction: Full matrix

### RIM-W

- RandGvecW: [RIM-W] Screened interaction RS components
- RandGvecW: [RIM-W] Screened interaction RS components
- rimw_type: [RIM-W] Screened interaction limit metal/semiconductor/graphene
- CUTGeo   [CUT] Coulomb Cutoff geometry: box/cylinder/sphere/ws/slab X/Y/Z/XY
- CUTBox   [CUT] [au] Box sides
- CUTRadius: [CUT] [au] Sphere/Cylinder radius
- CUTCylLen: [CUT] [au] Cylinder length
- CUTwsGvec: [CUT] WS cutoff: number of G to be modified
- CUTCol_test: [CUT] Perform a cutoff test in R-space
- ### Sxc
- GreenFTresh: [GW] [o/o] Treshold to define the new zoomed energy range
- QPExpand   [F GW] The QP corrections are expanded all over the BZ
- GreenF2QP  [F GW] Use real axis Green`s function to define the QP
- OnMassShell: [F GW] On mass shell approximation

### Real Time dynamics

- RTBands    [RT] Bands
- TwoAlpha   [RT] $C_{nk} \approx \alpha*\Gamma_{nk}^2 2_{alpha}$
- GrKind     [RT] G-ret kind: Lorentzian (QP)/ Hyperbolic QP_secant (HS)
- RADLifeTime: [RT] Radiative life-time (if negative Yambo sets it equal to Phase_LifeTime in NL)
- RADmagnific: [RT] Radiative life-time magnificationRAD_magnification
- PhLifeTime [RT] Constant Dephasing TimePhase_LifeTime
- DephTRange [RT] Time range in which Dephasing is applied
- DephEThresh [RT] Threshold on the energy difference between two states to dephase them

####  *_* Dynamics
- RTstep      [RT] Real Time step length
- NETime      [RT] Simulation Time
- dTupdateTimeSet: [RT] Time for manual deltaT update
- dTupdateTime: [RT] Initial Time for deltaT update (active only if non-zero)
- dTupdateJump: [RT] Time betweem two deltaT updates
- dTupdateTresh: [RT][o/o] Treshold of deltaT updates
- dT_MAX      [RT] Maximum value for the time-dependent dT
- dT_SET      [RT] Prefactor for manual dT update
- Integrator  [RT] Integrator. Use keywords space separated  ( "EULER/EXPn/INV" "SIMPLE/RK2/RK4/HEUN" "RWA")
- IO_times=(/CARR_RT_IO_t%INTERVAL_time_INPUT,Gless_RESTART_RT_IO_t%INTERVAL_time_INPUT,OUTPUT_RT_IO_t%INTERVAL_time_INPUT/)
- IOtime      [RT] Time between two consecutive I/O (CARRIERs - GF - OUTPUT)
- IOCachetime [RT] Time between two consecutive (caching - I/O) of OBSERVABLES
- RTehEny     [RT] Electron-hole energy range
####   *_flags_*
- DephCVonly    [RT] Dephase only in the CV channel
- RTskipImposeN [RT] Conservation of N, dN  imposed by hand on-the-fly
- RTEvalEnergy  [RT] Energy variation computed on the fly
- RTEvalEntropy [RT] Entropy variation computed on the fly
- SaveGhistory  [RT] Save the history of the green function
- *_updates_*
- RTUpdateSOC     [RT] Update the SOC interaction
- RTUpdateE     [RT] Update the Enery levels on-the-fly
- RTEqScatt     [RT] Include Gamma0f0 term in scattering
- RTImpForMet   [RT] Impose structure optimized for metals
- RTzeroTempRef [RT] Use zero temperature Fermi districution as reference
- RTskipPHabs   [RT] Skip e-p Lifetimes due to phonon absorption
  *_scattering_*
- LifeExtrapolation    [RT] Skipped Lifetimes are extrapolated
- LifeExtrapSteps   [RT] Step length between and inside two consecutive groups of lifetimes
- RelaxTimeApprox    [RT] Skipped Lifetimes are extrapolated
- RTAtemp   [RT] Temperatures for relaxation time approximation
- RTAchem   [RT] Chemical potentials for relaxation time approximation
- ScattTresh [RT] Treshold on the eh energy to compute the scatteringRT_scatt_tresh

####  *_EE scattering_*
- PlasmaPerc [RT] Plasma approx (0-100): % of eh pair consideredPLASMA_redux_percent,
- EERimPerc  [RT] EE Double Grid (0-100): % of the points used in EE scattDbGd_EE_percent,
- RTskipImposeE   [RT] Conservation of E (e-e channel) imposed by hand on-the-fly

### ELPH

- MemTresh   [RT] Treshold on the decay of the retarded GFNE_MEM_treshold,
- UseDebyeE   [RT] Use a single Debye energy for all phonon modes
- RT_T_evol   [RT] Use a complete Time evolution instead of the CCA
- InducedField: [RT] Include induced field in coupling and current
- VelGaugeCorr: [RT] Correct the non local term of the pseudo with the vector potential

#### OLD/EXPERIMNETAL

- RTAveDeph    [RT] Dephasing for all elements not included in RTDePhMatrix 
- LifeFitTemp: [RT] Fit on the fly  lifetimes ratio to a Fermi distribution

### NL

- NLBands      [NL] Bands range
- NLverbosity  [NL] Verbosity level (low | high)
- NLstep       [NL] Time step lengthRT_step,unit=Time_unit(1)
- NLtime       [NL] Simulation Time
- NLintegrator [NL] Integrator ("EULEREXP/RK2/RK4/RK2EXP/HEUN/INVINT/CRANKNIC")
- NLCorrelation: [NL] Correlation ("IPA/HARTREE/TDDFT/LRC/LRW/JGM/SEX")
- NLLrcAlpha   [NL] Long Range Correction
- NLDamping    [NL] Damping (or dephasing)
- NLEnRange    [NL] Energy range
- NLEnSteps    [NL] Energy steps
- UseDipoles: [NL] Use Covariant Dipoles (just for test purpose)
- FrSndOrd: [NL] Force second order in Covariant Dipoles
- EvalCurrent: [NL] Evaluate the current
- NLSampleWF  [NL] Sample the WFs (sampling NL order)n_order

### RT OR NL

if SC

- ExtF_Dir             [NL ExtF] Versor
- ExtF_Int             [NL ExtF] Intensity
- ExtF2_Dir             [NL ExtF] Versor
- ExtF2_Int             [NL ExtF] Intensity
- FrSndOrd: [NL] Force second order in Covariant Dipoles


SLKdim  [SLK] Matrix Dimension for scalapack

- Hamiltonian   [MAG] Hamiltonian kind [pauli,landau,all]MAG_hamiltonian_type
- B_Field       [MAG] Magnetic field modulus
- B_psi         [MAG] Magnetic field psi angle
- B_theta       [MAG] Magnetic field theta angle
- B_Gauge       [MAG] Gauge ("SYMM"etric, "X_ASYMM", "Y_ASYMM" or "Z_ASYMM"etric)
- PhaseTrick: [MAG] Phase trick for a better diagonalization
- B_radius      [MAG] Magnetic flux radius

#### *_Memory_*

- MEM_tresh     [MEMORY] Treshold on traced memory allocations/deallocations

### PHEL

- PHDbGdsize [PHEL] Size of subset of double grid k-points
- PHELQpts   [PHEL] Q-points considered
- PHELTrans  [PHEL] Energy window around W_ph to select transitions (units of GDamping)
- PHEL_QPH_En: [PHEL] Energy points to get the Quasi-Phonon solution (units of the bare PH energy)
- PH_SE_mode [PHEL] Self-Energy scattering mode ("bare-bare","dressed-bare","dressed-dressed")

- FxcMEStps: [TDDFT] [o/o] Memory energy steps
- FxcSVdig: [TDDFT] Keep SV that are within FxcSVdig digits
- FxcRetarded [TDDFT] Retarded TDDFT kernel

- BSehWind: [BSK] [o/o] E/h coupling pairs energy windowBS_eh_win,
- BSEQptR [BSK] Transferred momenta range
- BDmERef [BSS] Damping energy reference

### ACFDT

- AC_n_LAM [ACFDT] Coupling Costant GL-grid points
- AC_n_FR  [ACFDT] Integration frequency points
- AC_E_Rng [ACFDT] Imaginary axis 2nd & 3rd energy points

### DIPOLES

- DipBands       [DIP] Bands range for dipoles
- DipQpts        [DIP] Qpts range for dipoles
- DipoleEtresh   [DIP] Treshold in the definition of R=P
- DipApproach    [DIP] [G-space v/R-space x/Covariant/Shifted grids]
- DipComputed    [DIP] [default R P V; extra P2 Spin Orb]
- ShiftedPaths   [DIP] Shifted grids paths (separated by a space)grid_paths,
- DipPDirect [DIP] Directly compute $<v>$ also when using other approaches for dipoles
- DipBandsALL: [DIP] Compute all bands range, not only valence and conduction

### NL OR ELECTRIC

- EvPolarization: [DIP] Evaluate Polarization (require DipApproach=Covariant)
- FrPolPerdic: [DIP] Force periodicity of polarization respect to the external field

### RT

- SPINprojected [DIP] Project the spin dipoles in the c/v channels 

### ELPH

- GphBRnge  [ELPH] G[W] bands range
- ElPhModes [ELPH] Phonon modes
- FANdEtresh: [ELPH] Energy treshold for Fan denominator
- DWdEtresh [ELPH] Energy treshold for DW denominator
- GkkpDB    [ELPH] GKKP database (gkkp | gkkp_expanded | genFroh )
  
- ElPhHBRnge: [ELPH] Hamiltonian bands
- ElPhHKpt  [ELPH] Hamiltonian k-point
- REStresh  [ELPH] Residual treshold to report in output
- WRgFsq: [ELPH] Dump on file

### SC

- SCBands   [SC] Bands

- SCIter    [SC] SC Iterations
- SCRhoTresh: [SC] Rho convergence threshold
- SC_precondition: [SC] Kind of preconditionin: thomas-fermi, simple, none
- SCUpWIter [SC] Update W(q,G,G) every SCUpWIter iteractions
- Mean_Potential: [SC] Real-space Mean Potential
- SCnlMix: [SC] Use SC non-local mixing
- FrozeDensity: [NL] Do not update density (for testing purposes)
- SCneqKIND  [SC] Options are [contrained-occ/constrained-mu/matsubara]SC_neq_kind
- SCmu       [SC] Reference / holes / electrons chem potential
- SCcohIt    [SC] Impose off-diagonal rho in the initial basis set for N iterations
#### *_common with RT_*
- BandMix   [SC] Band mixingSC_band_mixing
- SCmixing  [SC] SC Cycle Mixing (< 1.)SC_cycle_mixing,
- SCdiag: [SC] Diagonal approximation for the self-energy(WF unchaged)
- SCEtresh  [SC] Energy convergence threshold for SC-GW

### SC OR RT OR QED

- COLLBands   [COLL] Bands for the collisions
- HXC_Potential [SC] SC HXC PotentialH_potential
- COLLCut    [SC,RT] Cutoff on the collisions
- OEPItSolver: [SC] Iterative solution instead of inversion of OEP
- OEPapprox: [SC] OEP approximation: n=none s=Slater k=KLI c=CED +w=Weighted