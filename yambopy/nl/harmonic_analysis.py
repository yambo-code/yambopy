# Copyright (c) 2023-2025, Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate linear response from real-time calculations (yambo_nl)
#
import numpy as np
from yambopy.units import ha2ev,fs2aut, SVCMm12VMm1,AU2VMm1
from yambopy.nl.external_efield import Divide_by_the_Field
from tqdm import tqdm
import scipy.linalg
import sys
import os

def Coefficients_Inversion(NW,NX,P,W,T_period,T_range,T_step,efield,INV_MODE):
    """
    Compute coefficients inversion using various inversion modes.
    see Sec. III in PRB 88, 235113 (2013) 

    Parameters:
        NW (int): Number of components (assumed NW = NX), order of the response functions.
        NX (int): Number of components.
        P (array): real-time polarizations.
        W (array): Frequency components, multiples of the laser frequency.
        T_period (float): Period of time.
        T_range (tuple): Start and end of time range.
        T_step (float): Time step.
        efield (dict): Contains external field information.
        INV_MODE (str): Inversion mode ('full', 'lstsq', 'svd').

    Returns:
        tuple: X_here (inverted coefficients), Sampling (array of time and values).
    """
    #
    # Here we use always NW=NX
    #
    M_size = 2*(NW-1) + 1  # Positive and negative components plut the zero
    nP_components = NX
    # 
    i_t_start = int(np.round(T_range[0]/T_step)) 
    i_deltaT  = int(np.round(T_period/T_step)/M_size)


# Memory alloction 
    M      = np.zeros((M_size, M_size), dtype=np.cdouble)
    X_here = np.zeros(nP_components, dtype=np.cdouble)


# Calculation of  T_i and P_i
    T_i = np.array([(i_t_start + i_deltaT * i) * T_step - efield["initial_time"] for i in range(M_size)])
    P_i = np.array([P[i_t_start + i_deltaT * i] for i in range(M_size)])
    Sampling = np.column_stack((T_i / fs2aut, P_i))

# Build the M matrix
    M[:, 0] = 1.0
    for i_n in range(1, nP_components):
        exp_neg = np.exp(-1j * W[i_n] * T_i, dtype=np.cdouble)
        exp_pos = np.exp(1j * W[i_n] * T_i, dtype=np.cdouble)
        M[:, i_n] = exp_neg
        M[:, i_n - 1 + NX] = exp_pos

# Invert M matrix
    INV_MODES = ['full', 'lstsq', 'svd']
    if INV_MODE not in INV_MODES:
        raise ValueError("Invalid inversion mode. Expected one of: %s" % INV_MODES)
  
    if INV_MODE=="full":
        try:
# Invert M matrix
            INV = np.linalg.inv(M)
        except:
            print("Singular matrix!!! standard inversion failed ")
            print("set inversion mode to LSTSQ")
            INV_MODE="lstsq"

    if INV_MODE=='lstsq':
# Least-squares
        I = np.eye(M_size,M_size)
        INV = np.linalg.lstsq(M, I, rcond=tol)[0]

    if INV_MODE=='svd':
# Truncated SVD
        INV = np.linalg.pinv(M,rcond=tol)

# Calculate X_here
    X_here=np.zeros(nP_components,dtype=np.cdouble)
    for i_n in range(nP_components):
        X_here[i_n]=X_here[i_n]+np.sum(INV[i_n,:]*P_i[:])

    return X_here,Sampling



def Harmonic_Analysis(nldb, X_order=4, T_range=[-1, -1],prn_Peff=False,INV_MODE="full",prn_Xhi=True):
    """
    Perform harmonic analysis on a dataset.

    Parameters:
        nldb: Input dataset with fields, polarizations, and time points.
        X_order (int): Maximum harmonic order.
        T_range (list): Time range for analysis.
        prn_Peff (bool): Print effective polarizations if True.
        INV_MODE (str): Inversion mode ('full', 'lstsq', 'svd').
        prn_Xhi (bool): Print susceptibilities if True.

    Returns:
        tuple: Frequencies and susceptibilities if prn_Xhi is False.
    """
        # Time series and step
    time = nldb.IO_TIME_points
    T_step = time[1] - time[0]
    efield = nldb.Efield[0]
    n_runs = len(nldb.Polarization)
    polarization = nldb.Polarization

    freqs = np.array([efield["freq_range"][0] for efield in nldb.Efield], dtype=np.double)

    print("\n* * * Harmonic analysis * * *\n")

    # Check for valid field type
    if efield["name"] not in {"SIN", "SOFTSIN", "ANTIRES"}:
        print("Harmonic analysis works only with SIN or SOFTSIN fields")
        sys.exit(0)

    # Check for single field
    if any(ef["name"] != "none" for ef in nldb.Efield_general[1:]):
        print("Harmonic analysis works only with a single field, please use sum_frequency.py functions")
        sys.exit(0)

    print(f"Number of runs: {n_runs}")

    W_step = min(freqs)
    max_W = max(freqs)
    
    print(f"Minimum frequency: {W_step * ha2ev:.3e} [eV]")
    print(f"Maximum frequency: {max_W * ha2ev:.3e} [eV]")

    T_period = 2.0 * np.pi / W_step
    print(f"Effective max time period: {T_period / fs2aut:.3f} [fs]")

    if T_range[0] <= 0.0:
        T_range[0] = time[-1] - T_period
    if T_range[1] <= 0.0:
        T_range[1] = time[-1]

    print(f"Time range: {T_range[0] / fs2aut:.3f} - {T_range[1] / fs2aut:.3f} [fs]")
    T_range_initial = np.copy(T_range)
        
    M_size = 2*X_order + 1  # Positive and negative components plut the zero
    X_effective       =np.zeros((X_order+1,n_runs,3),dtype=np.cdouble)
    Susceptibility    =np.zeros((X_order+1,n_runs,3),dtype=np.cdouble)
    Sampling          =np.zeros((M_size, 2,n_runs,3),dtype=np.double)    
    Harmonic_Frequency=np.zeros((X_order+1,n_runs),dtype=np.double)

    # Generate multiples of each frequency
    for i_order in range(X_order+1):
        Harmonic_Frequency[i_order,:]=i_order*freqs[:]
    
    loop_on_angles = nldb.n_angles != 0
    loop_on_frequencies = nldb.n_frequencies != 0

    if loop_on_angles:
        angles = np.linspace(0, 360, n_runs, endpoint=False)
        print("Loop on angles...")

    if loop_on_frequencies:
        print("Loop on frequencies...")
        
    # Find the Fourier coefficients by inversion
    for i_f in tqdm(range(n_runs)):
        T_period = 2.0 * np.pi / Harmonic_Frequency[1, i_f]
        T_range, out_of_bounds = update_T_range(T_period, T_range_initial, time)

        if out_of_bounds:
            print(f"WARNING! Time range out of bounds for frequency: {Harmonic_Frequency[1, i_f] * ha2ev:.3e} [eV]")

        for i_d in range(3):
            X_effective[:, i_f, i_d], Sampling[:, :, i_f, i_d] = Coefficients_Inversion(
                X_order + 1, X_order + 1, polarization[i_f][i_d, :],
                Harmonic_Frequency[:, i_f], T_period, T_range, T_step, efield, INV_MODE
            )
    # Calculate susceptibilities
    for i_order in range(X_order + 1):
        for i_f in range(n_runs):
            if i_order == 1:
                Susceptibility[i_order, i_f, 0] = 4.0 * np.pi * np.dot(efield['versor'], X_effective[i_order, i_f, :])
            else:
                Susceptibility[i_order, i_f, :] = X_effective[i_order, i_f, :]
            Susceptibility[i_order, i_f, :] *= Divide_by_the_Field(nldb.Efield[i_f], i_order)

    if nldb.calc != 'SAVE':
        prefix = f'-{nldb.calc}'
    else:
        prefix = ''


    # Reconstruct effective polarization
    if prn_Peff:
        print("Reconstruct effective polarizations...")
        Peff = np.zeros((n_runs, 3, len(time)), dtype=np.cdouble)
        for i_f in tqdm(range(n_runs)):
            for i_d in range(3):
                for i_order in range(X_order + 1):
                    freq_term = np.exp(-1j * i_order * freqs[i_f] * time)
                    Peff[i_f, i_d, :] += X_effective[i_order, i_f, i_d] * freq_term
                    Peff[i_f, i_d, :] += np.conj(X_effective[i_order, i_f, i_d]) * np.conj(freq_term)

        print("Print effective polarizations...")
        for i_f in tqdm(range(n_runs)):
            values = np.column_stack((time / fs2aut, Peff[i_f, 0, :].real, Peff[i_f, 1, :].real, Peff[i_f, 2, :].real))
            output_file = f'o{prefix}.YamboPy-pol_reconstructed_F{i_f + 1}'
            np.savetxt(output_file, values, header="[fs] Px Py Pz", delimiter=' ', footer="Reconstructed polarization")


    #Units of measure rescaling
    for i_order in range(X_order+1):
        Susceptibility[i_order,:,:]*=get_Unit_of_Measure(i_order)
    
    # Write final results
    if prn_Xhi:
        print("Write final results: xhi^1, xhi^2, xhi^3, etc...")
        for i_order in range(X_order + 1):
            output_file = f'o{prefix}.YamboPy-X_probe_order_{i_order}'
            header = "[eV] " + " ".join([f"X/Im(z){i_order} X/Re(z){i_order}" for _ in range(3)])
            values = np.column_stack((freqs * ha2ev, Susceptibility[i_order, :, 0].imag, Susceptibility[i_order, :, 0].real,
                                      Susceptibility[i_order, :, 1].imag, Susceptibility[i_order, :, 1].real,
                                      Susceptibility[i_order, :, 2].imag, Susceptibility[i_order, :, 2].real))
            np.savetxt(output_file, values, header=header, delimiter=' ', footer="Harmonic analysis results")
    else:
        return freqs, Susceptibility


def get_Unit_of_Measure(i_order):
    """
    Calculates the unit of measure based on the order provided.

    Parameters:
    i_order (int): The order that determines the calculation of the unit of measure.

    Returns:
    float: The calculated unit of measure.
    """
    # Define the constant ratio for conversion
    ratio = SVCMm12VMm1 / AU2VMm1

    # Calculate the unit of measure based on the order
    if i_order == 0:
        return ratio
    return np.power(ratio, i_order - 1, dtype=np.float64)

def update_T_range(t_period, t_range_initial, time):
    """
    Updates the time range for analysis based on the specified period.
    
    Parameters:
    t_period (float): The period for the analysis range.
    t_range_initial (list or tuple): Initial time range as [start, end].
    time (list or array): Array of time points for reference.

    Returns:
    tuple: Updated time range as [start, end], and a boolean indicating if the range is out of bounds.
    """
    # Initialize the analysis range
    t_range = list(t_range_initial)  # Ensure we work on a mutable copy
    t_range[1] = t_range[0] + t_period

    # Check if the range goes out of bounds and adjust if necessary
    t_range_out_of_bounds = False
    if t_range[1] > time[-1]:
        t_range[1] = time[-1]
        t_range[0] = t_range[1] - t_period
        t_range_out_of_bounds = True

    return t_range, t_range_out_of_bounds

