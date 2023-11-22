# Code from Myrta Gruning (2023)

import sys
import numpy as np
import matplotlib.pyplot as plt
from yambopy.units import fs2aut
from scipy import signal

def zeropadding_signal(f,a,Nval):
    step = f[1]-f[0]
    pf = np.zeros(2*Nval+len(f))
    pf[Nval-1:Nval+len(f)-1] = f[:]
    for n in range(Nval-1,-1,-1):
        pf[n] = pf[n+1] - step
    for n in range(Nval+1):
        pf[Nval+len(f)+n-1] = pf[Nval+len(f)+n-2] + step
#    pa = np.pad(a,(Nval,),'constant',constant_values=(a[0],a[-1])) #here I supposed that the first and last values are zero
    pa = np.pad(a,(Nval,),'constant',constant_values=(0,0)) #here I supposed that the first and last values are zero
    return pf,pa


def plot_signal(data=None,time=None,hdir=[1,0,0],padded=False,Npad=600,tstring='Signal',singlefig = True):
    #
    if data is None or time is None:
        print('Error you should provide an array with the data and time')
        sys.exit(0)
    
    hdir=hdir/np.linalg.norm(hdir)
    hdata=np.dot(data.T, hdir)

    if (singlefig):
        plt.show()
    plt.title(tstring)
    plt.xlabel('Time [fs]')
    if padded:
        nfreq,spadded=zeropadding_signal(time/fs2aut,hdata,Npad)
        plt.plot(nfreq,spadded)
    else:
        plt.plot(time/fs2aut,hdata)
    return

def get_psd(data=None, time=None,wind="blackman",hdir=[1,0,0],padded=False,Npad=600):
    #
    hdir=hdir/np.linalg.norm(hdir)
    hdata=np.dot(data.T, hdir)
    #
    if padded:
        nfreq,spadded=zeropadding_signal(time/fs2aut,hdata,Npad)
        f, Pxx_den = signal.periodogram(spadded,1.0/(nfreq[1]-nfreq[0]),scaling='spectrum',window=wind)
    else:
        f, Pxx_den = signal.periodogram(hdata,fs=fs2aut/(time[1]-time[0]),scaling='spectrum',window=wind)

# 1 hertz [Hz] = 4.13566553853599E-15 electron-volt [eV]

    return f*4.1356655385, Pxx_den/np.max(Pxx_den)

def plot_psd(f,psd,lfreq,tstring='',singlefig = True, lmax=41,ymin=10e-7, gap=None):
    l=np.arange(1,lmax)
    if (singlefig):
        plt.show()
        tlabel='Laser freq = '+str(lfreq)+'eV'
        ttitle = tstring
    else:
        tlabel = tstring
        ttitle = ''
    print(f)
    plt.semilogy(f/lfreq, psd,label=tlabel)
    plt.title(ttitle)
    plt.ylim([ymin,1])
    plt.xlim([0.1,lmax])
    plt.vlines(l,ymin,1,linestyles='dotted')
    plt.xticks(l)
    if gap is not None:
        plt.vlines(gap/lfreq,ymin,1,linestyles='dotted',colors='red')
    plt.xlabel('Harmonic number')
    plt.ylabel(r"PSD $\left[V^2\right]$")
    plt.legend()
    return

