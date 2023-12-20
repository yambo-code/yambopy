import numpy as np
import sys

def Fourier_Interpolation(ft, fw, time,freqs,mode="T2W"):
    #
    # I assume constant time-step and regular frequency grid
    # otherwise this subroutine needs to be generalized
    #
    t_step=time[1]-time[0]
    f_step=freqs[1]-freqs[0]
    #
    if mode.upper() == "T2W":
        fw[:,:]=0.0+0.0j
        for i_w in range(fw.shape[1]):
            for i_c in range(ft.shape[0]):
                    fw[i_c,i_w] = fw[i_c,i_w]+np.sum(ft[i_c,:]*np.exp(1j*freqs[i_w]*time[:]))*t_step
    elif mode.upper() == "W2T":
        ft[:,:]=0.0+0.0j
        for i_w,i_t in zip(range(ft.shape[1]),range(fw.shape[1])):
            for i_c in range(ft.shape[0]):
                ft[i_c,i_t] = ft[i_c,i_t]+fw[i_c,i_w]*np.exp(-1j*freqs[i_w]*time[i_t])*f_step


