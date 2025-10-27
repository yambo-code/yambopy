# Copyright (c) 2025, Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
# Generate external field for (yambo_rt/yambo_nl)
#
import math
import numpy as np
from yambopy.units import speed_of_light,speed_of_light_SI

def Efield_strength(Intensity, unit_system):
    """
    From Octopus ( http://www.tddft.org/programs/octopus )

    It is very common to describe the strength of a laser field by its intensity, 
    rather than using the electric field amplitude. In atomic units, the relationship 
    between instantaneous electric field and intensity is:

        I(t) = c / (8π) * E^2(t)

    It is common to read intensities in kW cm^-2. The dimensions of intensities are [W]/(L^2T), 
    where [W] are the dimensions of energy. The relevant conversion factors are:

        Hartree / (a_0^2 * atomic_time) = 6.4364086e+15 W / cm^2 = 6.4364086e+12 kW / cm^2

    In Yambo AU2KWCMm2 = 6.4364086e+12

    --------------------------------------------------------------------------------------------
    This simple function uses the formula valid in SI and CGS to extract 
    the field intensity given in atomic units.
    """

    # Work space
    SPEED = 0.0
    I = 0.0

    # From Boyd, "Nonlinear Optics", 3rd edition, page 602–603
    # Assuming n=1 (epsilon=1)

    if unit_system == "SI":
        # Convert input intensity to correct SI-based scale
        I = Intensity * AU2KWCMm2 / WMm22KWCMm2

        SPEED = speed_of_light_SI

        # I = 1 / (FREE_SPACE_PERM * SPEED_OF_LIGHT) * |E|^2 
        Efield_strength = math.sqrt(I / FREE_SPACE_PERM / SPEED) * VMm12AU

    elif unit_system == "CGS":
        # Convert input intensity to correct CGS-based scale
        I = Intensity * AU2KWCMm2 / WMm22KWCMm2 / WMm22ERGCMm2SECm1

        SPEED = speed_of_light_SI * 100.0  # cm/sec

        # I = SPEED_OF_LIGHT / (4π) * |E|^2 
        Efield_strength = math.sqrt(I * 4.0 * math.pi / SPEED) * SVCMm12VMm1 * VMm12AU

    elif unit_system == "AU":
        # In atomic units, use direct formula
        I = Intensity
        SPEED = SPEED_OF_LIGHT

        # I = SPEED_OF_LIGHT / (4π) * |E|^2 
        Efield_strength = math.sqrt(I * 4.0 * math.pi / SPEED)

    else:
        raise ValueError(f"Unknown unit system: {unit_system}")

    return Efield_strength


def theta_function(T, step, order):
    """
    Element-wise version of the Fortran theta_function.
    
    Parameters
    ----------
    T : float or np.ndarray
        Input variable (can be scalar or array).
    step : float
        Step size.
    order : int
        Order of the derivative (0, 1, or 2).
        
    Returns
    -------
    float or np.ndarray
        Value of the theta function (or its derivative) evaluated at T.
    """
    T = np.asarray(T, dtype=float)
    theta = np.zeros_like(T)

    # Case where |T/step| > 1.05
    mask_far = np.abs(T / step) > 1.05
    if order == 0:
        theta[mask_far & (T > 0)] = 1.0
    # Higher orders remain zero by default
    if np.all(mask_far):
        return theta

    # Determine step_case
    step_cases = np.full(T.shape, "", dtype="<U4")
    x = T / step
    step_cases[np.abs(x + 1.0) < 0.1] = "-1.0"
    step_cases[np.abs(x + 0.5) < 0.1] = "-0.5"
    step_cases[np.abs(x + 0.0) < 0.1] = " 0.0"
    step_cases[np.abs(x - 0.5) < 0.1] = "+0.5"
    step_cases[np.abs(x - 1.0) < 0.1] = "+1.0"

    # Apply the cases
    for case in ["-1.0", "-0.5", " 0.0", "+0.5", "+1.0"]:
        mask = step_cases == case
        if not np.any(mask):
            continue
        if case == "-1.0":
            if order == 2:
                theta[mask] = 1.0 / step**2
        elif case == "-0.5":
            if order == 0:
                theta[mask] = 0.25
            elif order == 1:
                theta[mask] = 0.50 / step
            elif order == 2:
                theta[mask] = 1.00 / step**2
        elif case == " 0.0":
            if order == 0:
                theta[mask] = 0.5
            elif order == 1:
                theta[mask] = 1.0 / step
            elif order == 2:
                theta[mask] = 0.0
        elif case == "+0.5":
            if order == 0:
                theta[mask] = 0.75
            elif order == 1:
                theta[mask] = 0.50 / step
            elif order == 2:
                theta[mask] = -1.00 / step**2
        elif case == "+1.0":
            if order == 0:
                theta[mask] = 1.0
            elif order == 2:
                theta[mask] = -1.00 / step**2

    return theta


def Compose_Field(a_pot,t_start,t_step,t_range,e_amp,field_fname):
    """
    Write an external field for yambo_rt/yambo_nl

    Parameters:
        a_pot: vector potential A''(t), A'(t), A(t)
        t_step: Time-step [au]
        t_start: Starting time [au]
        t_range: Time range [au]
        field_fname: External field file name.

    Returns:
        write external field on file
    """

    theta=theta_function(t_range-t_start,t_step,0)
    delta=theta_function(t_range-t_start,t_step,1)
    signf=theta_function(t_range-t_start,t_step,2)

    A_coeff      = speed_of_light*e_amp
    A_vecpot     = A_coeff*(a_pot[0]*theta)
    A_vecpot_vel = A_coeff*(a_pot[1]*theta +a_pot[0]*delta)
    A_vecpot_acc = A_coeff*(a_pot[2]*theta +a_pot[1]*delta+a_pot[0]*signf)
