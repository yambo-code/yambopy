import numpy as np
import matplotlib.pyplot as plt

"""
Script to draw band plots.

It is called by generate_bands.

Arguments are:
    data, plt_type, show
    
    - data: [prefix,nkpoints,nbands,kstps,gaps,eigen,(shifted_eigen)] 
    - plt_type:
        -- 'bands': standard band plot with "Fulvio"-style format and layout 
    - out_name: string to be attached to output file name
    - erange: energy window centered at Fermi level (single value in eV, None means all energies)
    - show: determine if plt.show() is called at the end
    - print_data: if True, print plotted data in text format
"""

def plot_driver(data,plt_type,out_name=None,erange=None,show=True,print_data=True):
    """
    This function is like the __init__ function of a class.
    
    This script may be transformed into a class when new plot types are added
    """
    # Scissored energies are present
    if len(data)==9: prefix,nkpoints,nbands,kstps,KPTs_labels,points,gaps,eigs_noscissor,eigen = data
    # Scissored energies are not present
    elif len(data)==8: prefix,nkpoints,nbands,kstps,KPTs_labels,points,gaps,eigen = data
    # Data missing
    else: raise IndexError('Incorrect number of elements in plot data.')

    if eigen.shape[1] != nbands:
        if eigen.shape[1] == 2*nbands: eigen = np.array( [ eigen[:,:nbands],  eigen[:,nbands:] ]  ) 
        else: raise ValueError("Mismatch between nbands and eigenvalue array shape")

    if plt_type=='bands': 
        data_to_plot = prefix,nkpoints,nbands,kstps,KPTs_labels,points,gaps,eigen
        electron_dispersion_plot(data_to_plot,out_name,erange,show)
    if print_data: print_out_files(kstps,eigen,prefix,out_name)

def print_out_files(kstps,eigen,prefix,out_nm):
    """
    Print data as *.dat file
    """

    # Output name
    if out_nm is not None: out_file_dat = "%s_%s.dat"%(prefix,out_nm)
    else:                  out_file_dat = "%s_bands.dat"%prefix

    l_magnetic = (len(eigen.shape)==3)

    # Create array to print
    if not l_magnetic:
        Nk,Nb = eigen.shape
        to_prnt = np.zeros((Nk,Nb+1))
        to_prnt[:,0] =kstps
        to_prnt[:,1:]=eigen

        # Save array
        np.savetxt(out_file_dat,to_prnt,fmt='%.6f')

    if l_magnetic:
        Nspin,Nk,Nb = eigen.shape
        to_prnt_up = np.zeros((Nk,Nb+1))
        to_prnt_up[:,0] =kstps
        to_prnt_up[:,1:]=eigen[0]
        to_prnt_dn = np.zeros((Nk,Nb+1))
        to_prnt_dn[:,0] =kstps
        to_prnt_dn[:,1:]=eigen[1]

        # Save arrays
        np.savetxt('SPIN1_'+out_file_dat,to_prnt_up,fmt='%.6f')
        np.savetxt('SPIN2_'+out_file_dat,to_prnt_dn,fmt='%.6f')

def electron_dispersion_plot(data,out_nm,erange,show):
    """
    Standard band plot with "Fulvio"-style format and layout
    """
    from math import floor,ceil
    from matplotlib import rc,rcParams

    # Fonts in standard TeX style
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    # Data
    prefix,nkpoints,nbands,kstps,KPTs_labels,KPTs,gaps,eigen = data

    # Output name
    if out_nm is not None: out_file_pdf = "%s_%s.pdf"%(prefix,out_nm)
    else:                  out_file_pdf = "%s_bands.pdf"%prefix

    l_magnetic = (len(eigen.shape)==3)

    # Initial preparations
    if gaps[0]<1.e-6: yref=0.
    else:             yref=min(gaps)
    if not l_magnetic: Nk,Nb = eigen.shape
    if l_magnetic:  Ns,Nk,Nb = eigen.shape
    ylims = np.array([yref-erange,yref+erange])/2.
    xlims = [kstps[0],kstps[-1]]
    for il in range(len(KPTs_labels)):
        if KPTs_labels[il]=='G': KPTs_labels[il]=r'$\Gamma$'
    ## linewidths
    band_linewidth  = 1.5
    frame_linewidth = 3./4.*band_linewidth
    faint_linewidth = band_linewidth/3.
    rcParams['axes.linewidth'] = frame_linewidth
    ## yticks
    l_yticks = True
    N_ticks = floor(ylims[1])-ceil(ylims[0])
    if erange is None or erange>=20.: l_yticks=False
    if erange<20. and erange>=10.:
        yticks = [ ceil(ylims[0]) + tick for tick in range(N_ticks) ]
    if erange<10. and erange>=2.:
        yticks = [0.]
        tick = yticks[0]
        while tick<ylims[1]:
            tick +=0.5
            yticks.append(tick)
        tick = yticks[0]
        while tick>ylims[0]:
            tick -=0.5
            yticks.insert(0,tick)
    if erange<2.:
        yticks = [0.]
        tick = yticks[0]
        while tick<ylims[1]:
            tick +=0.1
            yticks.append(tick)
        tick = yticks[0]
        while tick>ylims[0]:
            tick -=0.1
            yticks.insert(0,tick)
    if l_yticks: yticklabels = [str(ytick) for ytick in yticks]
    ## xticks
    xticks = [0.]+KPTs
    xticklabels = KPTs_labels

    # Start plot
    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_ylim(ylims)
    ax.set_ylabel(r'Energy (eV)',size=20)
    if l_yticks:
        ax.set_yticks(yticks[1:])
        ax.set_yticklabels(yticklabels[1:],size=20)
    ax.set_xlim(xlims)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels,size=20)
    ax.tick_params(direction='in',width=frame_linewidth,length=6,left=True,right=True)
    for point in KPTs[:-1]: ax.axvline(point,color='black',linewidth=frame_linewidth)
    if yref!=0.: ax.axhline(min(gaps),color='gray',linestyle='--',linewidth=faint_linewidth)
    ax.axhline(0,color='black',linewidth=faint_linewidth)

    # Draw plot
    if not l_magnetic:
        for ib in range(Nb): ax.plot(kstps,eigen[:,ib],ls='-',lw=band_linewidth,c='red')
    if l_magnetic:
        clrs = ['red','blue']
        for i_s in range(Ns):
            for ib in range(Nb): ax.plot(kstps,eigen[i_s][:,ib],ls='-',lw=band_linewidth,c=clrs[i_s])

    plt.savefig(out_file_pdf)
    if show: plt.show()
