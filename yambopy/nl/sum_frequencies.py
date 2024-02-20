# Copyright (c) 2023, Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate linear response from real-time calculations (yambo_nl)
#
import numpy as np
from yambopy.units import ha2ev,fs2aut, SVCMm12VMm1,AU2VMm1
from yambopy.nl.external_efield import Divide_by_the_Field
from yambopy.nl.harmonic_analysis import update_T_range
import scipy.linalg
import sys
import os

#
# Polarization coefficient inversion see Sec. III in PRB 88, 235113 (2013) 
#
#  NW          order of the response functions 
#  NX          numer of coefficents required
#  P           real-time polarization 
#  W           multiples of the laser frequency
#  T_prediod   shorted cicle period
#  X           coefficents of the response functions X1,X2,X3...
#
def SF_Coefficents_Inversion(NW,NX,P,W1,W2,T_period,T_range,T_step,efield,tol,INVMODE="full"):
    #
    # Here we use always NW=NX
    #
    M_size = (2*(NW-1) + 1)**2  # Positive and negative components plut the zero
    # 
    i_t_start = int(np.round(T_range[0]/T_step)) 
    i_deltaT  = int(np.round(T_period/T_step)/M_size)


# Memory alloction 
    M      = np.zeros((M_size, M_size), dtype=np.cdouble)
    P_i    = np.zeros(M_size, dtype=np.double)
    T_i    = np.zeros(M_size, dtype=np.double)

# Calculation of  T_i and P_i
    for i_t in range(M_size):
        T_i[i_t] = (i_t_start + i_deltaT * i_t)*T_step - efield["initial_time"]
        P_i[i_t] = P[i_t_start + i_deltaT * i_t]

# Build the M matrix
    C = np.zeros((M_size, M_size), dtype=np.int8)
    for i_t in range(M_size):
        i_c = 0
        for i_n in range(-NX+1, NX):
            for i_n2 in range(-NX+1, NX):
                M[i_t, i_c]          = np.exp(-1j * (i_n*W1+i_n2*W2) * T_i[i_t],dtype=np.cdouble)
                C[i_n+NX-1,i_n2+NX-1] = i_c
                i_c+=1

# Multiple possibilities to calculate the inversion
    INV_MODES = ['full', 'lstsq', 'svd']
        if INVMODE not in INV_MODES:
            raise ValueError("Invalid inversion mode. Expected one of: %s" % INV_MODES)

    if INVMODE=="full":
        try:
# Invert M matrix
            INV = np.zeros((M_size, M_size), dtype=np.cdouble)
            INV = np.linalg.inv(M)
    except:
            print("Singular matrix!!! standard inversion failed ")
            print("set inversion mode to LSTSQ")
            INVMODE="lstsq"

    if INVMODE=='lstsq':
# Least-squares
        I = np.eye(M_size,M_size)
        INV = np.linalg.lstsq(M, I, rcond=tol)[0]

    if INVMODE=='svd'
# Truncated SVD
        INV = np.zeros((M_size, M_size), dtype=np.cdouble)
        INV = np.linalg.pinv(M,rcond=tol)

# Calculate X_here
    X_here=np.zeros((M_size, M_size),dtype=np.cdouble)
    for i_n in range(-NX+1, NX):
        for i_n2 in range(-NX+1, NX):
            for i_t in range(M_size):
                X_here[i_n+NX-1,i_n2+NX-1]=X_here[i_n+NX-1,i_n2+NX-1]+INV[C[i_n+NX-1,i_n2+NX-1],i_t]*P_i[i_t]

    return X_here



def SF_Harmonic_Analysis(nldb, tol=1e-7, X_order=4, T_range=[-1, -1],prn_Peff=False):
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
        print("Harmonic analysis works only with SIN or SOFTSIN fields")
        sys.exit(0)

    l_test_one_field=False
    if(nldb.Efield_general[1]["name"] == "SIN" or nldb.Efield_general[1]["name"] == "SOFTSIN"):
        # frequency of the second and third laser, respectively)
        pump_probe=nldb.Efield_general[1]["freq_range"][0] 
        print("Frequency of the second field : "+str(pump_probe*ha2ev)+" [eV] \b")
    elif(nldb.Efield_general[1]["name"] == "none"):
        print("Only one field present, please use standard harmonic_analysis.py for SHG,THG, etc..!")
        print(" * * * Test mode with a single field * * * ")
        print(" * * * Frequency of the second field assumed to be zero * * *")
        l_test_one_field=True
        pump_probe=0.0
    else:
        print("Fields different from SIN/SOFTSIN are not supported ! ")
        sys.exit(0)
    
    if(nldb.Efield_general[2]["name"] != "none"):
        print("Three fields not supported yet ! ")
        sys.exit(0)


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

    M_size = (2*X_order + 1)**2
    X_effective       =np.zeros((M_size,M_size,n_frequencies,3),dtype=np.cdouble)
    Susceptibility    =np.zeros((M_size,M_size,n_frequencies,3),dtype=np.cdouble)
    

    # Find the Fourier coefficients by inversion
    for i_f in range(n_frequencies):
        #
        # T_period change with the laser frequency 
        #
        #T_period=2.0*np.pi/Harmonic_Frequency[1,i_f]  #Not sure about this
        T_range,T_range_out_of_bounds=update_T_range(T_period,T_range_initial,time)
        #
        if T_range_out_of_bounds:
            print("WARNING! Time range out of bounds for frequency :",Harmonic_Frequency[1,i_f]*ha2ev,"[eV]")
        #
        for i_d in range(3):
            X_effective[:,:,i_f,i_d]=SF_Coefficents_Inversion(X_order+1, X_order+1, polarization[i_f][i_d,:],freqs[i_f],pump_probe,T_period,T_range,T_step,efield,tol)

    # Calculate Susceptibilities from X_effective
    for i_order in range(-X_order,X_order+1):
        for i_order2 in range(-X_order,X_order+1):
            for i_f in range(n_frequencies):
                Susceptibility[i_order+X_order,i_order2+X_order,i_f,:]=X_effective[i_order+X_order,i_order2+X_order,i_f,:]
            
                Susceptibility[i_order+X_order,i_order2+X_order,i_f,:]*=Divide_by_the_Field(nldb.Efield[i_f],abs(i_order)+abs(i_order2))

    if(prn_Peff):
        # Print time dependent polarization
        P=np.zeros((n_frequencies,3,len(time)),dtype=np.cdouble)
        for i_f in range(n_frequencies):
            for i_d in range(3):
                for i_order in range(-X_order,X_order+1):
                    for i_order2 in range(-X_order,X_order+1):
                        P[i_f,i_d,:]+=X_effective[i_order+X_order,i_order2+X_order,i_f,i_d]*np.exp(-1j * (i_order*freqs[i_f]+i_order2*pump_probe) * time[:])
        header2="[fs]            "
        header2+="Px     "
        header2+="Py     "
        header2+="Pz     "
        footer2='Time dependent polarization reproduced from Fourier coefficients'
        for i_f in range(n_frequencies):
            values2=np.c_[time.real/fs2aut]
            values2=values2=np.append(values2,np.c_[P[i_f,0,:].real],axis=1)
            values2=values2=np.append(values2,np.c_[P[i_f,1,:].real],axis=1)
            values2=values2=np.append(values2,np.c_[P[i_f,2,:].real],axis=1)
            output_file2='o.YamboPy-pol_reconstructed_F'+str(i_f+1)
            np.savetxt(output_file2,values2,header=header2,delimiter=' ',footer=footer2)

    # Print the result
    for i_order in range(-X_order,X_order+1):
        for i_order2 in range(-X_order,X_order+1):

            if i_order==0 and i_order2==0: 
                Unit_of_Measure = SVCMm12VMm1/AU2VMm1
            else:
                Unit_of_Measure = np.power(SVCMm12VMm1/AU2VMm1,abs(i_order)+abs(i_order2)-1,dtype=np.double)
        
            Susceptibility[i_order+X_order,i_order2+X_order,:,:]=Susceptibility[i_order+X_order,i_order2+X_order,:,:]*Unit_of_Measure

            output_file='o.YamboPy-SF_probe_order_'+str(i_order)+'_'+str(i_order2)
            if i_order == 0 or (i_order == 1 and i_order2 == 0) or (i_order == 0 and i_order2 == 1):
                header="E [eV]            X/Im(x)            X/Re(x)            X/Im(y)            X/Re(y)            X/Im(z)            X/Re(z)"
            else:
                header="[eV]            "
                header+="X/Im[cm/stV]^%d     X/Re[cm/stV]^%d     " % (abs(i_order)+abs(i_order2)-1,abs(i_order)+abs(i_order2)-1)
                header+="X/Im[cm/stV]^%d     X/Re[cm/stV]^%d     " % (abs(i_order)+abs(i_order2)-1,abs(i_order)+abs(i_order2)-1)
                header+="X/Im[cm/stV]^%d     X/Re[cm/stV]^%d     " % (abs(i_order)+abs(i_order2)-1,abs(i_order)+abs(i_order2)-1)

            values=np.c_[freqs*ha2ev]
            values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order,:,0].imag],axis=1)
            values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order,:,0].real],axis=1)
            values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order,:,1].imag],axis=1)
            values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order,:,1].real],axis=1)
            values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order,:,2].imag],axis=1)
            values=np.append(values,np.c_[Susceptibility[i_order+X_order,i_order2+X_order,:,2].real],axis=1)

            footer='Non-linear response analysis performed using YamboPy'
            np.savetxt(output_file,values,header=header,delimiter=' ',footer=footer)
