# Code from Myrta Gruning (2023)

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def zeropadding_signal(f,a,Nval):
    step = f[1]-f[0]
    pf = np.zeros(2*Nval+len(f))
    pf[Nval-1:Nval+len(f)-1] = f[:]
    for n in range(Nval-1,-1,-1):
        pf[n] = pf[n+1] - step
    for n in range(Nval+1):
        pf[Nval+len(f)+n-1] = pf[Nval+len(f)+n-2] + step
    pa = np.pad(a,(Nval,),'constant',constant_values=(a[0],a[-1])) #here I supposed that the first and last values are zero
    return pf,pa


def plot_signal(fname=None,array=None,time=None,idir=1,padded=False,Npad=600,tstring='Signal',singlefig = True):
    #
    if fname is not None:
        print('Ciao ciao')
        data=np.genfromtxt(fname,comments="#")
    elif array is not None and time is not None:
        data=array
    else: 
        print('Error you should provide a file or an array with the data')
        sys.exit(0)

    if (singlefig):
        plt.show()
    plt.title(tstring)
    plt.xlabel('Time [fs]')
    if padded:
        nfreq,spadded=zeropadding_signal(data[:,0],data[:,idir],Npad)
        plt.plot(nfreq,spadded)
    else:
        plt.plot(data[:,0],data[:,idir])
    return

def get_psd(fname,wind="blackman",idir=1,padded=False,Npad=600):
    #
    data=np.genfromtxt(fname,comments="#")
    if padded:
        nfreq,spadded=zeropadding_signal(data[:,0],data[:,idir],Npad)
        f, Pxx_den = signal.periodogram(spadded,1.0/(nfreq[1]-nfreq[0]),scaling='spectrum',window=wind)
    else:
        f, Pxx_den = signal.periodogram(data[:,idir],1.0/(data[1,0]-data[0,0]),scaling='spectrum',window=wind)
    return f*4.1356655385, Pxx_den/np.max(Pxx_den)

def plot_psd(f,psd,lfreq,tstring='',singlefig = True):
    l=np.arange(1,31)
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
    plt.ylim([10e-14,1])
    plt.xlim([0.1,30])
    plt.vlines(l,10e-14,1,linestyles='dotted')
    plt.vlines(2.26/lfreq,10e-14,1,linestyles='dotted',colors='red')
    plt.xlabel('Harmonic number')
    plt.ylabel(r"PSD $\left[V^2\right]$")
    plt.legend()
    return

