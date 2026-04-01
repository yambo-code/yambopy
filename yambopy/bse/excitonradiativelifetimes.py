# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: MRF
#
# This file is part of the yambopy project
#
import numpy as np
from yambopy import YamboLatticeDB,YamboExcitonDB,YamboDipolesDB
from yambopy.units import kb,autime2s,ha2ev,m_e,speed_of_light
import os
from yambopy.tools.citations import citation

def get_exciton_dipole(state, blongdir, ylat, ydip, yexc, _cache={}):
    """
    This function computes the dipole D of exciton state `state`
    at Q=0 starting from the transition-space expression:

    D_{state} = \sum_{cvk} A^{state}_{cvk} Efield \cdot r_cvk 

    :: A      --> BSE eigenvectors
    :: Efield --> electric field
    :: r      --> dipole matrix elements

    This is useful to change Efield direction without rerunning Yambo

    NB: in Yambo, q=0 is actually q=q0_def_norm=1e-5, so if you are checking for 
        consistency with the residuals in BS_diago db you have to multiply 
        by q0_def_norm

    Input:
    * blongdir:   electric field polarization direction
    * lattice_path, dipoles_path, bse_path : paths to required databases

    Output:
    * Values of |D|^2 in bohr^2
    """

    blongdir = np.asarray(blongdir, dtype=float)
    key = (id(ydip.dipoles), id(yexc.table), tuple(blongdir))

    # cache dipoles_proj and kcv_idx for same inputs
    if key not in _cache:
        field_dir = blongdir / np.linalg.norm(blongdir)
        dipoles_proj = np.einsum('x,kxcv->kcv', field_dir, ydip.dipoles)

        table_kcv = yexc.table[:, [0, 2, 1]]
        kcv_idx = tuple((table_kcv[:, :3] - 1).T)

        trans_dipoles = dipoles_proj[kcv_idx]
        _cache[key] = trans_dipoles
    else:
        trans_dipoles = _cache[key]

    dip_exc = np.dot(yexc.eigenvectors[state], trans_dipoles)
    dip_exc_squared = np.abs(dip_exc)**2

    return dip_exc_squared / ylat.nkpoints


@citation("H.-Y. Chen et al. Phys. Rev. B 100, 075135 (2019)")
def get_radiative_lifetime_3D_iso(T,state,ylat,ydip,yexc,Meff,eps,shiftE=0):

    """
        Function to compute the radiative lifetime tau_S(T) of a single exciton state for *isotropic* 3D materials.

        Input
            * T:        temperature in K
            * state:    index of the exciton state. Python convention (indices starting from 0)
            * ylat:     lattice database, YamboLatticeDB
            * ydip:     dipoles database, YamboDipolesDB
            * yexc:     BSE database, YamboExcitonDB
            * Meff:     exciton effective mass in electron rest mass (m_e) units. 
            * eps:      material's relative optical dielectric constant
            * shiftE:   rigid shift of exciton energies. Default: shift=0 eV

        Output
            * Radiative lifetime of exciton state in seconds
    """

    if not (isinstance(Meff,float) and isinstance(eps,float)):
        raise ValueError("Meff and eps must be floats")

    Omega = ylat.lat_vol  #Bohr**3

    muS2 = get_exciton_dipole(state,[1,1,1],ylat,ydip,yexc) # Bohr**2
    ES = np.real(yexc.eigenvalues[state])+shiftE  # eV

    Meff_eV = Meff*m_e   # eV
    
    T_factor = ( ES**2/(2*Meff_eV*kb*T) )**(3/2)
    dipole_factor = 8/3*np.sqrt(np.pi*eps)*muS2/Omega
    prefactor = 4*np.pi/autime2s
    
    return 1/(prefactor*dipole_factor*T_factor)

@citation("H.-Y. Chen et al. Phys. Rev. B 100, 075135 (2019)")
def get_radiative_lifetime_3D_aniso(T,state,ylat,ydip,yexc,Meff,eps,shiftE=0):
    """
        Function to compute the radiative lifetime tau_S(T) of a single exciton state for *anisotropic* (uniaxial) 3D materials.

        Input
            * T:        temperature in K
            * state:    index of the exciton state. Python convention (indices starting from 0)
            * ylat:     lattice database, YamboLatticeDB
            * ydip:     dipoles database, YamboDipolesDB
            * yexc:     BSE database, YamboExcitonDB
            * Meff:     2 element list: exciton effective masses in electron rest mass (m_e) units in-plane, out-of-plane: [Meff_xy, Meff_z]
            * eps:      2 element list: material's relative optical dielectric constants in-plane, out-of-plane: [eps_xy, eps_z]
            * shiftE:   rigid shift of exciton energies. Default: shift=0 eV

        Output
            * Radiative lifetime of exciton state in seconds
    """

    if not (isinstance(Meff,list) and len(Meff) == 2):
        raise ValueError("Meff must be lists of 2 floats: [Meff_xy, Meff_z]")
    if not (isinstance(eps,list) and len(eps) == 2):
        raise ValueError("eps must be lists of 2 floats: [eps_xy, eps_z]")

    Omega = ylat.lat_vol  #Bohr**3

    muS2_xy = get_exciton_dipole(state,[1,1,0],ylat,ydip,yexc) # Bohr**2
    muS2_z = get_exciton_dipole(state,[0,0,1],ylat,ydip,yexc)  # Bohr**2
    ES = np.real(yexc.eigenvalues[state])+shiftE  # eV

    Meff_xy_eV = Meff[0]*m_e  # eV
    Meff_z_eV = Meff[1]*m_e   # eV

    eps_xy = eps[0]
    eps_z = eps[1]

    T_factor = ( ES**2/(2*Meff_xy_eV**(2/3)*Meff_z_eV**(1/3)*kb*T) )**(3/2)
    dipole_factor = np.sqrt(np.pi*eps_xy)/Omega * ((2*eps_z/(3*eps_xy)+2)*muS2_xy+8/3*muS2_z)
    prefactor = 4*np.pi/autime2s
    
    return 1/(prefactor*dipole_factor*T_factor)  # seconds

@citation("H.-Y. Chen et al. Phys. Rev. B 100, 075135 (2019)")
def get_radiative_lifetime_2D(T,state,ylat,ydip,yexc,Meff,eps=1,dip_dir=[1,1,0],shiftE=0):
    """
        Function to compute the radiative lifetime tau_S(T) of a single exciton state for 2D materials.

        Input
            * T:        temperature in K
            * state:    index of the exciton state. Python convention (indices starting from 0)
            * ylat:     lattice database, YamboLatticeDB
            * ydip:     dipoles database, YamboDipolesDB
            * yexc:     BSE database, YamboExcitonDB
            * Meff:     exciton effective mass in electron rest mass (m_e) units. 
            * eps:      environment dielectric constant. Default: eps=1, vacuum.
            * dip_dir:  orientation of the dipole moments. Default: [1,1,0] in the xy plane.
            * shiftE:   rigid shift of exciton energies. Default: shift=0 eV

        Output
            * Radiative lifetime of exciton state in seconds
    """

    if not (isinstance(Meff,float) and isinstance(eps,float)):
        raise ValueError("Meff and eps must be floats")

    lat = ylat.lat
    A = np.cross(lat[0],lat[1])[2] # Bohr**2

    muS2 = get_exciton_dipole(state,dip_dir,ylat,ydip,yexc) # Bohr**2
    ES = np.real(yexc.eigenvalues[state])+shiftE # eV

    Meff_eV = Meff*m_e   # eV
    
    T_factor = 4/3*(ES**2/(2*Meff_eV*kb*T))
    gamma0_factor = (ES/ha2ev)*muS2/(eps*A)
    prefactor = 4*np.pi/speed_of_light/autime2s

    return 1/(prefactor*gamma0_factor*T_factor)   # seconds

@citation("H.-Y. Chen et al. Phys. Rev. B 100, 075135 (2019)")
def get_radiative_lifetime_1D(T,state,ylat,ydip,yexc,Meff,eps=1,dip_dir=[0,0,1],shiftE=0):
    """
        Function to compute the radiative lifetime tau_S(T) of a single exciton state for 1D materials.

        Input
            * T:        temperature in K
            * state:    index of the exciton state. Python convention (indices starting from 0)
            * ylat:     lattice database, YamboLatticeDB
            * ydip:     dipoles database, YamboDipolesDB
            * yexc:     BSE database, YamboExcitonDB
            * Meff:     exciton effective mass in electron rest mass (m_e) units. 
            * eps:      environment dielectric constant. Default: eps=1, vacuum.
            * dip_dir:  orientation of the dipole moments. Default: [0,0,1] along z.
            * shiftE:   rigid shift of exciton energies. Default: shift=0 eV

        Output
            * Radiative lifetime of exciton state in seconds
    """

    if not (isinstance(Meff,float) and isinstance(eps,float)):
        raise ValueError("Meff and eps must be floats")

    lat = ylat.lat
    Lz = lat[2,2]

    muS2 = get_exciton_dipole(state,dip_dir,ylat,ydip,yexc) # Bohr**2
    ES = np.real(yexc.eigenvalues[state])+shiftE  # eV

    Meff_eV = Meff*m_e   # eV
    
    T_factor = 4/3*(ES**2/(2*Meff_eV*kb*T))
    gamma0_factor = (ES/ha2ev)**2*muS2/(eps*Lz)
    prefactor = 4*np.pi/speed_of_light**2/autime2s
    
    return 1/(prefactor*gamma0_factor*T_factor)   # seconds

@citation("H.-Y. Chen et al. Phys. Rev. B 100, 075135 (2019)")
def get_radiative_lifetime_0D(state,ylat,ydip,yexc,eps=1,shiftE=0):
    """
        Function to compute the radiative lifetime tau_S(T) of a single exciton state for 0D materials.

        Input
            * state:    index of the exciton state. Python convention (indices starting from 0)
            * ylat:     lattice database, YamboLatticeDB
            * ydip:     dipoles database, YamboDipolesDB
            * yexc:     BSE database, YamboExcitonDB
            * eps:      environment dielectric constant. Default: eps=1, vacuum
            * shiftE:   rigid shift of exciton energies. Default: shift=0 eV.

        Output
            * Radiative lifetime of exciton state in seconds
    """

    if not isinstance(eps,float):
        raise ValueError("eps must be float")

    muS2 = get_exciton_dipole(state,[1,1,1],ylat,ydip,yexc) # Bohr**2
    ES = np.real(yexc.eigenvalues[state])+shiftE  # eV
    
    gamma0_factor = np.sqrt(eps)*(ES/ha2ev)**3*muS2
    prefactor = 4*np.pi/(3*speed_of_light**3)/autime2s
    
    return 1/(prefactor*gamma0_factor)   # seconds



def average_lifetime(Trange,states,ylat,ydip,yexc,dimension,Meff=None,eps=None,dip_dir=None,shiftE=0):
    """
        Function to compute the material's exciton radiative lifetime thermally averaged on the exciton states `states`.

        Input
            * Trange:       temperature range, in K
            * states:       list of indices of the exciton states to consider. Python convention (indices starting from 0.
            * ylat:         lattice database, YamboLatticeDB
            * ydip:         dipoles database, YamboDipolesDB
            * yexc:         BSE database, YamboExcitonDB
            * dimension:    dimensionality of the system. It can be:
                            - '3D' or '3D_iso' for isotropic 3D bulk systems
                            - '3D_aniso' for anisotropic (uniaxial) 3D bulks
                            - '2D', '1D', '0D' for lower-dimensional systems
            * Meff:         exciton effective mass. Needed by 3D, 2D, 1D calculations.
            * eps:          environment dielectric constant.
            * dip_dir:      dipole orientation. Needed by 2D and 1D calculations.
            * shiftE:   rigid shift of exciton energies. Default: shift=0 eV.

        Output
            * Array of averaged radiative lifetime with length as Trange in input
    """
    
    ES = np.real(yexc.eigenvalues)+shiftE
    DeltaE = ES-ES[states[0]]

    if dimension in ['3D','3D_iso']:
        gammas = np.array([np.exp(-DeltaE[state]/(kb*Trange))/get_radiative_lifetime_3D_iso(Trange,state,ylat,ydip,yexc,Meff,eps,shiftE) for state in states])
    elif dimension=='3D_aniso':
        gammas = np.array([np.exp(-DeltaE[state]/(kb*Trange))/get_radiative_lifetime_3D_aniso(Trange,state,ylat,ydip,yexc,Meff,eps,shiftE) for state in states])
    elif dimension=='2D':
        gammas = np.array([np.exp(-DeltaE[state]/(kb*Trange))/get_radiative_lifetime_2D(Trange,state,ylat,ydip,yexc,Meff,eps,dip_dir,shiftE) for state in states])
    elif dimension=='1D':
        gammas = np.array([np.exp(-DeltaE[state]/(kb*Trange))/get_radiative_lifetime_1D(Trange,state,ylat,ydip,yexc,Meff,eps,dip_dir,shiftE) for state in states])
    elif dimension=='0D':
        gammas = np.array([np.exp(-DeltaE[state]/(kb*Trange))/get_radiative_lifetime_0D(state,ylat,ydip,yexc,eps,shiftE) for state in states])
    else:
        raise ValueError("dimension must be one of the following: '3D','3D_iso','3D_aniso', '2D', '1D', '0D'")

    gamma_sum = np.sum(gammas,axis=0)

    normalize = np.array([np.exp(-DeltaE[state]/(kb*Trange)) for state in states])
    normalize_sum = np.sum(normalize,axis=0)

    return 1/(gamma_sum/normalize_sum)  # seconds

