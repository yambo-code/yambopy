# Copyright (c) 2023, Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate linear response from real-time calculations (yambo_nl)
#
import numpy as np
from yambopy.units import ha2ev,fs2aut, SVCMm12VMm1,AU2VMm1
from yambopy.nl.external_efield import Divide_by_the_Field
import scipy.linalg
import sys
import os
from fractions import Fraction
from sklearn.decomposition import TruncatedSVD
from scipy.linalg import svd
from numpy import zeros
from numpy import diag

#
# Polarization coefficient inversion see Sec. III in PRB 88, 235113 (2013) 
#
#  NW          order of the response functions 
#  NX          number of coefficents required
#  P           real-time polarization 
#  W           multiples of the laser frequency
#  T_period    shorted cicle period
#  X           coefficents of the response functions X1,X2,X3...
#
def Coefficents_Inversion(NW,NX,P,W,T_period,T_range,T_step,efield):
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
    P_i    = np.zeros(M_size, dtype=np.double)
    T_i    = np.zeros(M_size, dtype=np.double)
    X_here = np.zeros(nP_components, dtype=np.cdouble)

# Calculation of  T_i and P_i
    for i_t in range(M_size):
        T_i[i_t] = (i_t_start + i_deltaT * i_t)*T_step - efield["initial_time"]
        P_i[i_t] = P[i_t_start + i_deltaT * i_t]

# Build the M matrix
    for i_t in range(M_size):
        M[i_t, 0] = 1.0

    for i_t in range(M_size):
        for i_n in range(1, nP_components):
            M[i_t, i_n]          = np.exp(-1j * W[i_n] * T_i[i_t],dtype=np.cdouble)
            M[i_t, i_n - 1 + NX] = np.exp( 1j * W[i_n] * T_i[i_t],dtype=np.cdouble)

# Truncated SVD
    #svd=TruncatedSVD(n_components=M_size)
    #svd.fit(M)
    #result=svd.inverse_transform(M)
    
    #U, s, VT = svd(M)
    #Sigma = zeros((M.shape[0], M.shape[1]))
    #Sigma[:M.shape[1], :M.shape[1]] = diag(s)
    #B2 = U.dot(Sigma.dot(VT))
    #print(B2-M)
    #Sigma = zeros((M.shape[0], M.shape[1]))
    #Sigma[:M.shape[0], :M.shape[0]] = diag(s)
    #n_elements = M_size
    #Sigma = Sigma[:, :n_elements]
    #VT = VT[:n_elements, :]
    #B = U.dot(Sigma.dot(VT))

# Invert M matrix
    INV = scipy.linalg.inv(M)
    #INV2 = scipy.linalg.inv(B)
    #INV_test= np.transpose(VT).dot(Sigma.dot(np.transpose(U)))
    #print(INV2-INV_test)
    #print(B-M)

# Calculate X_here
    X_here=np.zeros(nP_components,dtype=np.cdouble)
    for i_n in range(nP_components):
        for i_t in range(M_size):
           X_here[i_n]=X_here[i_n]+INV[i_n,i_t]*P_i[i_t] 

    return X_here



def Harmonic_Analysis(nldb, Order, T_range=[-1, -1]):
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

    freqs=np.zeros(n_frequencies,dtype=np.double)

    if efield["name"] != "SIN" and efield["name"] != "SOFTSIN" and efield["name"] != "ANTIRES":
        print("Harmonic analysis works only with SIN or SOFTSIN fields")
        sys.exit(0)

    print("\n* * * Harmonic analysis * * *\n")

    print("Number of frequencies : %d " % n_frequencies)
    # Smaller frequency
    W_step=sys.float_info.max
    max_W =sys.float_info.min
    for count, efield in enumerate(nldb.Efield):
        freqs[count]=efield["freq_range"][0]
        if efield["freq_range"][0]<W_step:
            W_step=efield["freq_range"][0]
        if efield["freq_range"][0]>max_W:
            max_W=efield["freq_range"][0]
    print("Minimum frequency : ",str(W_step*ha2ev)," [eV] ")
    print("Maximum frequency : ",str(max_W*ha2ev)," [eV] ")

    
    # Period of the incoming laser
    T_period=2.0*np.pi/W_step
    print("Effective max time period ",str(T_period/fs2aut)+" [fs] ")

    if T_range[0] <= 0.0:
        T_range[0]=time[-1]-T_period
    if T_range[1] <= 0.0:
        T_range[1]=time[-1]
    
    T_range_initial=np.copy(T_range)

    print("Time range : ",str(T_range[0]/fs2aut),'-',str(T_range[1]/fs2aut),'[fs]')

    f=np.vectorize(int)
    output_order=np.zeros(n_frequencies,dtype=np.int16)
    for i_f in range(n_frequencies):
        output_order[i_f]=f(np.around((freqs[i_f]*ha2ev+2*1.16)/0.036250))
    X_order=max(output_order)+Order
    #print(output_order)
    #print(X_order)

    X_effective       =np.zeros((X_order+1,n_frequencies,3),dtype=np.cdouble)
    Susceptibility    =np.zeros((n_frequencies,3),dtype=np.cdouble)
    Harmonic_Frequency=np.zeros((X_order+1),dtype=np.double)

    # Generate multiples of each frequency
    for i_order in range(X_order+1):
        Harmonic_Frequency[i_order]=i_order*freqs[0]

    #print(Harmonic_Frequency[output_order]*ha2ev)

    # Find the Fourier coefficients by inversion
    for i_f in range(n_frequencies):
        #
        # T_period change with the laser frequency 
        #
        T_period=2.0*np.pi/Harmonic_Frequency[1]
        #print(T_period/fs2aut)
        T_range,T_range_out_of_bounds=update_T_range(T_period,T_range_initial,time)
        #
        if T_range_out_of_bounds:
            print("WARNING! Time range out of bounds for frequency :",Harmonic_Frequency[1]*ha2ev,"[eV]")
        #
        for i_d in range(3):
            X_effective[:,i_f,i_d]=Coefficents_Inversion(X_order+1, X_order+1, polarization[i_f][i_d,:],Harmonic_Frequency,T_period,T_range,T_step,efield)
        Susceptibility[i_f,:]=X_effective[output_order[i_f],i_f,:]
        Susceptibility[i_f,:]*=Divide_by_the_Field(nldb.Efield[i_f],3)   # here X_order=3 since we have three fields always

        
    # Print time depenent polarization
    P=np.zeros((n_frequencies,3,len(time)),dtype=np.cdouble)
    for i_f in range(n_frequencies):
        for i_d in range(3):
            P[i_f,i_d,:]=X_effective[0,i_f,i_d]
            for i_order in range(1,X_order+1):
                P[i_f,i_d,:]=P[i_f,i_d,:]+X_effective[i_order,i_f,i_d].real*np.cos(Harmonic_Frequency[i_order] * time[:],dtype=np.cdouble)+X_effective[i_order,i_f,i_d].imag*np.sin(Harmonic_Frequency[i_order] * time[:],dtype=np.cdouble)
    for i_f in range(n_frequencies):
        values2=np.c_[time.real/fs2aut]
        values2=np.append(values2,np.c_[P[i_f,0,:].real],axis=1)
        values2=np.append(values2,np.c_[P[i_f,1,:].real],axis=1)
        values2=np.append(values2,np.c_[P[i_f,2,:].real],axis=1)
        output_file2='o.Pol_'+str(i_f+1)
        header2="[fs]            "
        header2+="Px     "
        header2+="Py     "
        header2+="Pz     "
        footer2='Time dependent polarization reproduced from Fourier coefficients'
        np.savetxt(output_file2,values2,header=header2,delimiter=' ',footer=footer2)


    # Print the result for CARS
    Unit_of_Measure = np.power(SVCMm12VMm1/AU2VMm1,2,dtype=np.double)
    Susceptibility[:,:]=Susceptibility[:,:]*Unit_of_Measure
    output_file='o.YamboPy-X_probe_CARS'
    header="[eV]            "
    header+="X/Im[cm/stV]^2     X/Re[cm/stV]^2     "
    header+="X/Im[cm/stV]^2     X/Re[cm/stV]^2     "
    header+="X/Im[cm/stV]^2     X/Re[cm/stV]^2     "
    values=np.c_[freqs*ha2ev]
    values=np.append(values,np.c_[Susceptibility[:,0].imag],axis=1)
    values=np.append(values,np.c_[Susceptibility[:,0].real],axis=1)
    values=np.append(values,np.c_[Susceptibility[:,1].imag],axis=1)
    values=np.append(values,np.c_[Susceptibility[:,1].real],axis=1)
    values=np.append(values,np.c_[Susceptibility[:,2].imag],axis=1)
    values=np.append(values,np.c_[Susceptibility[:,2].real],axis=1)
    footer='Non-linear response analysis performed using YamboPy'
    np.savetxt(output_file,values,header=header,delimiter=' ',footer=footer)
    

def update_T_range(T_period,T_range_initial,time):
        #
        # Define the time range where analysis is performed
        #
        T_range=T_range_initial
        T_range[1]=T_range[0]+T_period
        #
        # If the range where I perform the analysis is out of bounds
        # I redefine it
        #
        T_range_out_of_bounds=False
        #
        if T_range[1] > time[-1]:
            T_range[1]= time[-1]
            T_range[0]= T_range[1]-T_period
            T_range_out_of_bounds=True

        return T_range,T_range_out_of_bounds

    

