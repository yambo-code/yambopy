# Copyright (c) 2023, Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate linear response from real-time calculations (yambo_nl)
#
import numpy as np
from yambopy.nl.fft_interp import *
from yambopy.nl.external_efield import *
from yambopy.units import ha2ev,fs2aut

from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import os

def sci_format(x,lim):
    return '{:.1e}'.format(x)

def Plot_Pol_or_Curr(time=None, pol=None, curr=None, xlim=None,save_file=None):
    if not isinstance(pol, np.ndarray) and isinstance(pol, np.ndarray):
        print("Polarzation or Current not present")
        return
    if not isinstance(time, np.ndarray):
        print("Time series not present")
        return

    char_size=14


    fig=plt.figure(figsize=(10,10))
    fig.subplots_adjust(hspace=.05)

    major_formatter = FuncFormatter(sci_format)

    if xlim == None:
        xlim=[0.0,time[-1]/fs2aut]

    if isinstance(pol, np.ndarray):
        pj='P'
        arr=pot
        fig.suptitle(' Real-time polarization in the three cartesian directions ', fontsize=char_size)
    else:
        pj='J'
        arr=curr
        fig.suptitle(' Real-time current in the three cartesian directions ', fontsize=char_size)

    ax=plt.subplot(3,1,1)
    plt.tick_params(axis='both', which='major', labelsize=char_size)
    plt.plot(time[:]/fs2aut,arr[0,:],color="blue", linewidth=1.5, linestyle="-",label=pj+"$_x$(t)")
    ax.yaxis.set_major_formatter(major_formatter)
    ax.set_xticks([])
    ax.set_xlim(xlim)
    ax.legend(fontsize=char_size)

    ###########6
    ax=plt.subplot(3,1,2)
    plt.tick_params(axis='both', which='major', labelsize=char_size)
    plt.plot(time[:]/fs2aut,arr[1,:],color="blue", linewidth=1.5, linestyle="-",label=pj+"$_y$(t)")
    ax.yaxis.set_major_formatter(major_formatter)
    ax.set_xticks([])
    ax.set_xlim(xlim)
    ax.legend(fontsize=char_size)
    
    ax=plt.subplot(3,1,3)
    plt.tick_params(axis='both', which='major', labelsize=char_size)
    plt.plot(time[:]/fs2aut,arr[2,:],color="blue", linewidth=1.5, linestyle="-",label=pj+"$_z$(t)")
    ax.set_xlabel('[fs]',size=14)
    ax.set_xlim(xlim)
    ax.legend(fontsize=char_size)
    if save_file is not None:
        plt.savefig(save_file)

    plt.show()

def Linear_Response(time, pol, efield, pol_ref=None, e_range=[0.0, 20.0], n_freqs=200, 
        output_file="o.YamboPy-eps_along_E",plot=False,plot_file=None):
    if efield["name"] != "DELTA":
        print("Linear response implemented only for Delta function external fields ")
        sys.exit(0)

    # Pump and probe minus the pump if present
    if pol_ref is not None:
        pol=pol-pol_ref

    freqs=np.linspace(e_range[0],e_range[1],n_freqs)/ha2ev
    pol_w=np.zeros((pol.shape[0],n_freqs),dtype=complex)

    efield_w=get_Efield_w(freqs,efield)

    Fourier_Interpolation(pol,pol_w,time,freqs,mode="T2W")
    #
    # Polarization along the field direction
    #
    pol_w_along_E=np.zeros(n_freqs,dtype=complex)
    for i_d in range(3):
        pol_w_along_E[:]+=pol_w[i_d,:]*efield["versor"][i_d]
    #
    # EPS = 1+4.0*pi*Xhi(omega)
    #
    eps=np.zeros(n_freqs,dtype=complex)
    eps=1.0+4.0*np.pi*pol_w_along_E/efield_w

    header="E [eV]      Im/eps      Re/eps"
    footer='Linear response analysis performed using YamboPy'
    np.savetxt(output_file,np.c_[freqs*ha2ev,eps.imag,eps.real],header=header,delimiter=' ',footer=footer)

    fig, axs = plt.subplots(2)
    fig.suptitle('Dielectric constant along the field direction')
    axs[0].set_xlim(e_range)
    axs[0].plot(freqs*ha2ev,eps.real,label='Real part ')
    axs[0].legend()
    axs[0].set_xticks([])
    axs[1].set_xlim(e_range)
    axs[1].plot(freqs*ha2ev,eps.imag,label='Imag part ')
    axs[1].legend()
    axs[1].set_xlabel('eV')
    if plot:
        plt.show()
    if plot_file is not None:
        plt.savefig(plot_file)
