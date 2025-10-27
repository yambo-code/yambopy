# Copyright (c) 2025, Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
# Generate external field for (yambo_rt/yambo_nl)
#
import numpy as np
from yambopy.units import speed_of_light

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


def Compose_Field(a_pot,t_start,t_step,t_range,field_fname):
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

    A_coeff      = speed_of_light*
    A_vecpot     = a_pot[0]*theta
    A_vecpot_vel = a_pot[1]*theta +a_pot[0]*delta 
    A_vecpot_acc = a_pot[2]*theta +a_pot[1]*delta+a_pot[0]*signf
