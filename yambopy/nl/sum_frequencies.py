# Copyright (c) 2023, Mike Nico Pionteck and Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate linear response from real-time calculations (yambo_nl)
#
import numpy as np
import math
from yambopy.units import ha2ev,fs2aut, SVCMm12VMm1,AU2VMm1
from yambopy.nl.external_efield import Divide_by_the_Field
from yambopy.nl.harmonic_analysis import update_T_range
from scipy.optimize import least_squares
from tqdm import tqdm
from scipy.ndimage import uniform_filter1d
import itertools

import scipy.linalg
import sys
import os

#
# Polarization coefficient inversion see Sec. III in PRB 88, 235113 (2013) 
#
#  N_samp      number of sampling points 
#  NX+1        numer of coefficents required
#  P           real-time polarization 
#  W           multiples of the laser frequency
#  T_prediod   shorted cicle period
#  X           coefficents of the response functions X1,X2,X3...
#
def SF_Coefficents_Inversion(N_samp,NX,NX2,P,W1,W2,T_range,T_step,efield,tol,INV_MODE,SAMP_MOD,INV0=None):
    #
    M_size = (2*NX+1)*(2*NX2+1)  # Positive and negative components plus the zero
    #
    if N_samp<=M_size: 
        raise ValueError(" Too few sampling points please increase it ")

    i_t_start = int(np.round( T_range[0] / T_step)) 
    i_deltaT  = int(np.round((T_range[1]-T_range[0])/T_step)/N_samp)

# Memory alloction 
    M        = np.zeros((N_samp, M_size), dtype=np.cdouble)
    P_i      = np.zeros(N_samp, dtype=np.double)
    T_i      = np.zeros(N_samp, dtype=np.double)
    Sampling = np.zeros((N_samp,2), dtype=np.double)

# Calculation of  T_i and P_i
    SAMP_MODES = {'linear', 'log', 'random'}
    if SAMP_MOD not in SAMP_MODES:  
        raise ValueError(f"Invalid sampling mode. Expected one of: {SAMP_MODES}")
    #
    if SAMP_MOD=='linear':
        T_i = (i_t_start + i_deltaT * range(N_samp))*T_step - efield["initial_time"]
        P_i = P[i_t_start + i_deltaT * range(N_samp)]
    elif SAMP_MOD=='log':
        T_i = np.geomspace(i_t_start * T_step, T_range[1], N_samp, endpoint=False)
        P_i = [P[int(np.round(t / T_step))] for t in T_i]
    elif SAMP_MOD=='random':
        T_i = np.random.uniform(i_t_start * T_step, T_range[1], N_samp)
        P_i = [P[int(np.round(t / T_step))] for t in T_i]
    
    Sampling[:,0]=T_i/fs2aut
    Sampling[:,1]=P_i

# Build the M matrix
    C = np.zeros((2*NX+1, 2*NX2+1), dtype=np.int8)
    for i_t in range(N_samp):
        for i_c,(i_n,i_n2) in enumerate(itertools.product(range(-NX, NX+1),range(-NX2, NX2+1))):
            M[i_t, i_c]          = np.exp(-1j * (i_n*W1+i_n2*W2) * T_i[i_t],dtype=np.cdouble)
            C[i_n+NX,i_n2+NX2] = i_c

# Multiple possibilities to calculate the inversion
    INV_MODES = {'full', 'lstsq', 'svd','lstsq_init'}
    if INV_MODE not in INV_MODES: 
        raise ValueError(f"Invalid inversion mode. Expected one of:  {INV_MODES} ")

    if INV_MODE=="full":
# Invert M matrix
        if N_samp != M_size: 
            raise TypeError("Only square matrix can be used with full inversion")
        INV = np.zeros((M_size, M_size), dtype=np.cdouble)
        try:
            INV = np.linalg.inv(M)
        except:
            print("Singular matrix!!! standard inversion failed ")
            print("set inversion mode to LSTSQ")
            INV_MODE="lstsq"
    if INV_MODE=='lstsq_init':
        def residuals_func(x):
            x_cmplx=x[0:int(x.size/2)] + 1j * x[int(x.size/2):x.size]
            return np.linalg.norm(np.dot(M, x_cmplx) - P_i)

        # This function works only with real values
        # I convert the complex x0 in to real
        if(INV0 is None):
            x0_cmplx = np.linalg.lstsq(M, P_i, rcond=tol)[0]
        else:
            x0_cmplx = INV0
        x0 = np.concatenate((x0_cmplx.real, x0_cmplx.imag))
        res = least_squares(residuals_func, x0, ftol=1e-11,gtol=1e-11,xtol=1e-11,verbose=1,x_scale='jac')
        INV = res.x[0:int(res.x.size/2)] + 1j * res.x[int(res.x.size/2):res.x.size]

    if INV_MODE=='lstsq':
# Least-squares
        INV = np.linalg.lstsq(M, P_i, rcond=tol)[0]

    if INV_MODE=='svd':
# Truncated SVD
        INV = np.zeros((N_samp, M_size), dtype=np.cdouble)
        INV = np.linalg.pinv(M,rcond=tol)

# Calculate X_here
    X_here=np.zeros((2*NX+1, 2*NX2+1),dtype=np.cdouble)
    for i_n,i_n2 in itertools.product(range(-NX, NX+1),range(-NX2, NX2+1)):
        i_c=C[i_n+NX,i_n2+NX2]
        if INV_MODE=='lstsq' or INV_MODE=='lstsq_init':
            X_here[i_n+NX,i_n2+NX2]=INV[i_c]
        else:
            X_here[i_n+NX,i_n2+NX2]=X_here[i_n+NX,i_n2+NX2]+np.sum(INV[i_c,:]*P_i[:])

    return X_here,Sampling,INV


def SF_Harmonic_Analysis(nldb, tol=1e-10, X_order=4, X_order2=None, T_range=[-1, -1], N_samp=-1,prn_Peff=False,prn_Xhi=True,INV_MODE='svd',SAMP_MOD='log'):
    # Time series 
    time  =nldb.IO_TIME_points
    # Time step of the simulation
    T_step=nldb.IO_TIME_points[1]-nldb.IO_TIME_points[0]
    # External field of the first run
    efield=nldb.Efield[0]
    # Numer of exteanl laser frequencies
    n_frequencies=len(nldb.Polarization)
    # Array of polarizations for each laser frequency
    polarization=nldb.Polarization

    print("\n* * * Sum/difference frequency generation: harmonic analysis * * *\n")

    freqs=np.zeros(n_frequencies,dtype=np.double)

    if efield["name"] != "SIN" and efield["name"] != "SOFTSIN" and efield["name"] != "ANTIRES":
        raise ValueError("Harmonic analysis works only with SIN or SOFTSIN fields")

    l_test_one_field=False

    if(X_order2==None): X_order2=X_order

    if(nldb.Efield_general[1]["name"] == "SIN" or nldb.Efield_general[1]["name"] == "SOFTSIN"):
        # frequency of the second and third laser, respectively)
        pump_freq=nldb.Efield_general[1]["freq_range"][0] 
        print("Frequency of the second field : "+str(pump_freq*ha2ev)+" [eV] \b")
    elif(nldb.Efield_general[1]["name"] == "none"):
        print("Only one field present, please use standard harmonic_analysis.py for SHG,THG, etc..!")
        print("* * * Test mode with a single field * * * ")
        print(" * * * Frequency of the second field assumed to be zero * * *")
        print(" * * * Only results for SHG,THG,... are correct * * * ")
        l_test_one_field=True
        pump_freq=0.0
    else:
        raise ValueError("Fields different from SIN/SOFTSIN are not supported ! ")
    
    if(nldb.Efield_general[2]["name"] != "none"):  raise ValueError("Three fields not supported yet ! ")

    print("Number of frequencies : %d " % n_frequencies)
    # Smaller frequency
    W_step=sys.float_info.max
    max_W =sys.float_info.min

    for count, efield in enumerate(nldb.Efield):
        freqs[count]=efield["freq_range"][0]
        
        if efield["freq_range"][0]<W_step: W_step=efield["freq_range"][0]
        if efield["freq_range"][0]>max_W:   max_W=efield["freq_range"][0]

    print("Minimum frequency : ",str(W_step*ha2ev)," [eV] ")
    print("Maximum frequency : ",str(max_W*ha2ev)," [eV] ")
    
    # Period of the incoming laser
    T_period=2.0*np.pi/W_step
    print("Effective max time period for field1 ",str(T_period/fs2aut)+" [fs] ")

    if T_range[0] <= 0.0: T_range[0]=2.0/nldb.NL_damping*6.0

    if T_range[1] <= 0.0: T_range[1]=time[-1]
    
    T_range_initial=np.copy(T_range)

    print("Initial time range : ",str(T_range[0]/fs2aut),'-',str(T_range[1]/fs2aut)," [fs] ")
    print("Pump frequency : ",str(pump_freq*ha2ev),' [eV] ')

    M_size = (2*X_order + 1)*(2*X_order2+1)

    if N_samp==-1: N_samp = M_size*2

    print(" Number of coefficents : "+str(M_size))
    print(" Number of sampling points : "+str(N_samp))

    X_effective       =np.zeros((2*X_order+1,2*X_order2+1,n_frequencies,3),dtype=np.cdouble)
    Sampling          =np.zeros((N_samp,2,n_frequencies,3),dtype=np.double)
    Susceptibility    =np.zeros((2*X_order+1,2*X_order2+1,n_frequencies,3),dtype=np.cdouble)
    INV0              =np.zeros((M_size,n_frequencies,3),dtype=np.cdouble)

    
    print("Loop in frequecies...")
    # Find the Fourier coefficients by inversion
#    old_tol=tol
    for i_f in tqdm(range(n_frequencies)):
#        These commented lines increase tol if two frequencies are degenerate
#        tol=old_tol
#        if np.isclose(pump_freq,freqs[i_f],atol=1e-7):
#            print(" WARNING: frequency "+str(i_f+1)+" = "+str(freqs[i_f]*ha2ev)+ " very close to the pump one: inversion tolerance reduced ")
#            tol=tol*100.0
        for i_d in range(3):
            X_effective[:,:,i_f,i_d],Sampling[:,:,i_f,i_d],INV0[:,i_f,i_d]=SF_Coefficents_Inversion(N_samp, X_order, X_order2, polarization[i_f][i_d,:],freqs[i_f],pump_freq,T_range,T_step,efield,tol,INV_MODE,SAMP_MOD)
        
        
# check non-converged points and degneracies and fix them
    spike_correction=False

    if(spike_correction):
        #
        # Response function to check for spike
        i_d=1
        i_order=1
        i_order2=1
        #
        iX=[i_order+X_order,i_order2+X_order2]
    # Calculate the moving average with a window
#       signal= abs(X_effective[i_order+X_order,i_order2+X_order,:,i_d])
        signal_im= abs(X_effective[i_order+X_order,i_order2+X_order2,:,i_d].imag)
        signal_re= abs(X_effective[i_order+X_order,i_order2+X_order2,:,i_d].real)

        window_size = 5
   #     smooth_signal = uniform_filter1d(signal, size=window_size)
        smooth_signal_im = uniform_filter1d(signal_im, size=window_size)
        smooth_signal_re = uniform_filter1d(signal_re, size=window_size)
    # Identify spikes relative to the local average
        threshold_local = 0.3  # defines how much a value can deviate from the local average
#        spike_indices_local = np.where(np.abs(signal - smooth_signal)/smooth_signal > threshold_local)[0]
        spike_indices_local_im = np.where(np.abs(signal_im - smooth_signal_im)/smooth_signal_im > threshold_local)[0]
        spike_indices_local_re = np.where(np.abs(signal_re - smooth_signal_re)/smooth_signal_re > threshold_local)[0]
        spike_indices_local = np.unique(np.concatenate((spike_indices_local_im, spike_indices_local_re)))

        print("Spike indices: ",spike_indices_local)
        for i_f in tqdm(spike_indices_local):
           if(i_f==0 and i_f in spike_indices_local):
               INV0[:,i_f,i_d]=INV0[:,i_f+1,i_d]
           elif(i_f==n_frequencies-1 and i_f in spike_indices_local):
               INV0[:,i_f,i_d]=INV0[:,i_f-1,i_d]
           else:
               INV0[:,i_f,i_d]=(INV0[:,i_f+1,i_d]+INV0[:,i_f-1,i_d])/2.0
               X_effective[:,:,i_f,i_d],Sampling[:,:,i_f,i_d],INV0[:,i_f,i_d]=SF_Coefficents_Inversion(N_samp, X_order, X_order2, polarization[i_f][i_d,:],freqs[i_f],pump_freq,T_range,T_step,efield,tol,INV_MODE="lstsq_init",SAMP_MOD=SAMP_MOD,INV0=INV0[:,i_f,i_d])

    print("Calculate susceptibility ")
    for i_order,i_order2 in itertools.product(range(-X_order,X_order+1),range(-X_order2,X_order2+1)):
        Susceptibility[i_order+X_order,i_order2+X_order2,:,:]=X_effective[i_order+X_order,i_order2+X_order2,:,:]
        if l_test_one_field:
            Susceptibility[i_order+X_order,i_order2+X_order2,:,:]*=Divide_by_the_Field(nldb.Efield[0],abs(i_order))
        else:
            D2=1.0
            if i_order!=0:
                D2*=Divide_by_the_Field(nldb.Efield[0],abs(i_order))
            if i_order2!=0:
                D2*=Divide_by_the_Field(nldb.Efield2[0],abs(i_order2))
            Susceptibility[i_order+X_order,i_order2+X_order2,:,:]*=D2

    if nldb.calc!='SAVE':
        prefix='-'+nldb.calc
    else:
        prefix=''

    if(prn_Peff):
        print("Reconstruct effective polarizations ...")        
        # Print time dependent polarization
        P=np.zeros((n_frequencies,3,len(time)),dtype=np.cdouble)
        for i_f,i_d in tqdm(itertools.product(range(n_frequencies),range(3))):
            for i_order,i_order2 in itertools.product(range(-X_order,X_order+1),range(-X_order2,X_order2+1)):
                P[i_f,i_d,:]+=X_effective[i_order+X_order,i_order2+X_order2,i_f,i_d]*np.exp(-1j * (i_order*freqs[i_f]+i_order2*pump_freq) * time[:])

        header2="[fs]            "
        header2+="Px     "
        header2+="Py     "
        header2+="Pz     "
        footer2='Time dependent polarization reproduced from Fourier coefficients'
        print("Write reconstructed polarizations ...")        
        for i_f in tqdm(range(n_frequencies)):
            values=np.c_[time.real/fs2aut]
            values=np.append(values,np.c_[P[i_f,0,:].real],axis=1)
            values=np.append(values,np.c_[P[i_f,1,:].real],axis=1)
            values=np.append(values,np.c_[P[i_f,2,:].real],axis=1)
            output_file2='o'+prefix+'.YamboPy-pol_reconstructed_F'+str(i_f+1)
            np.savetxt(output_file2,values,header=header2,delimiter=' ',footer=footer2)

        # Print Sampling point
        footer2='Sampled polarization'
        print("Write sampling ...")        
        for i_f in tqdm(range(n_frequencies)):
            values=np.c_[Sampling[:,0,i_f,0]]
            values=np.append(values,np.c_[Sampling[:,1,i_f,0]],axis=1)
            values=np.append(values,np.c_[Sampling[:,1,i_f,1]],axis=1)
            values=np.append(values,np.c_[Sampling[:,1,i_f,2]],axis=1)
            output_file3='o'+prefix+'.YamboPy-sampling_F'+str(i_f+1)
            np.savetxt(output_file3,values,header=header2,delimiter=' ',footer=footer2)

        print("Print general error in P(t) reconstruction ")
        footer2='Error in reconstructed polarization'
        header2="[eV]            "
        header2+="err[Px]     "
        header2+="err[Py]     "
        header2+="err[Pz]     "
        i_t_start = int(np.round(T_range[0]/T_step)) 
        values=np.zeros((n_frequencies,4),dtype=np.double)
        N=len(P[i_f,i_d,:])-i_t_start
        print("Write error ...")        
        for i_f in tqdm(range(n_frequencies)):
            values[i_f,0]=freqs[i_f]*ha2ev
            for i_d in range(3):
                values[i_f,i_d+1]=np.sqrt(np.sum((P[i_f,i_d,i_t_start:].real-polarization[i_f][i_d,i_t_start:]))**2)/N
        output_file4='o'+prefix+'.YamboPy-errP'
        np.savetxt(output_file4,values,header=header2,delimiter=' ',footer=footer2)
                

    # Print the result
    print("Write susceptibilities ...")        
    for i_order,i_order2 in itertools.product(range(-X_order,X_order+1),range(-X_order2,X_order2+1)):
        if i_order==0 and i_order2==0: 
            Unit_of_Measure = SVCMm12VMm1/AU2VMm1
        else:
            Unit_of_Measure = np.power(SVCMm12VMm1/AU2VMm1,abs(i_order)+abs(i_order2)-1,dtype=np.double)
            Susceptibility[i_order+X_order,i_order2+X_order2,:,:]=Susceptibility[i_order+X_order,i_order2+X_order2,:,:]*Unit_of_Measure
        output_file='o'+prefix+'.YamboPy-SF_probe_order_'+str(i_order)+'_'+str(i_order2)
        if i_order == 0 or (i_order == 1 and i_order2 == 0) or (i_order == 0 and i_order2 == 1):
            header="E [eV]            X/Im(x)            X/Re(x)            X/Im(y)            X/Re(y)            X/Im(z)            X/Re(z)"
        else:
            header="[eV]            "
            header+="X/Im[cm/stV]^%d     X/Re[cm/stV]^%d     " % (abs(i_order)+abs(i_order2)-1,abs(i_order)+abs(i_order2)-1)
            header+="X/Im[cm/stV]^%d     X/Re[cm/stV]^%d     " % (abs(i_order)+abs(i_order2)-1,abs(i_order)+abs(i_order2)-1)
            header+="X/Im[cm/stV]^%d     X/Re[cm/stV]^%d     " % (abs(i_order)+abs(i_order2)-1,abs(i_order)+abs(i_order2)-1)

        values=np.c_[freqs*ha2ev]
        values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order2,:,0].imag],axis=1)
        values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order2,:,0].real],axis=1)
        values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order2,:,1].imag],axis=1)
        values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order2,:,1].real],axis=1)
        values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order2,:,2].imag],axis=1)
        values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order2,:,2].real],axis=1)

        footer='Non-linear response analysis performed using YamboPy'
        if prn_Xhi:  
            np.savetxt(output_file,values,header=header,delimiter=' ',footer=footer)

    return Susceptibility,freqs


def update_T_range(T_range_initial,pump_freq, probe_freq):
    dec=1  # use only the first decimal in eV
    #
    a = round(pump_freq*ha2ev*10**dec)
    b = round(probe_freq*ha2ev*10**dec)
    r = a*b
    c = a*10**dec
    d = b*10**dec
    T_range=np.copy(T_range_initial)
    T_test=lcm(c,d)/r*ha2ev*2.0*np.pi+T_range[0]
    if T_test<T_range[1]:
        if round(a/b,3)==round(pump_freq/probe_freq,3):
             T_range[1]=T_test
        else:
            print("False multiple ",str(probe_freq*ha2ev))
# check for false multiply
    return T_range

def lcm(a,b):
  n = a
  m = b
  if (n < m):
    i = m
    m = n
    n = i
  p = n
  while p != 0:
    p = m%n
    m = n
    n = p
  return a*b/m
