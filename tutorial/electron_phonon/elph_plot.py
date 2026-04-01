"""
Tutorial for YamboElectronPhononDB and LetzElphElectronPhononDB.

Electron-phonon matrix element reading and plotting

EDIT the path below to point to the yambo SAVE folder.
"""
save_path='.'
from yambopy import YamboLatticeDB,YamboElectronPhononDB,LetzElphElectronPhononDB
from yambopy.units import ha2ev,bohr2ang
from yambopy.plot.plotting import shifted_grids_2D,BZ_hexagon
import numpy as np
import matplotlib.pyplot as plt

def get_plottable(elph,b0,b1,units='eV'):
    """
    Prepare the quantity to be plotted, which is:

    :: sqrt(  1/4 sum_{n1,n2} |g_{n1n2}^M_def (Q_def,k)|^2 )

    - M_def (mode) and Q_def (Q) already selected in input elph
    - sum factor is 1/4 because we assume only 2 consecutive bands, either b0=0,b1=1 or b0=2,b1=3
    - expression is converted from Rydberg to meV
    """
    if units=='meV':   dim = 1.
    elif units=='eV':  dim = 1000.
    elif units=='Ha':  dim = ha2ev*1000.
    elif units=='Ry':  dim = ha2ev/2.*1000.
    else: raise ValueError("elph units not recognized")

    elph2plot = np.sum( elph[:,b0:b1+1,b0:b1+1], axis=(1,2) )
    return np.sqrt(elph2plot/4.)*dim

def plot_2D_elph(pts,data,rlat=None,plt_cbar=False,**kwargs):
    """
    This is a custom version of the `plot_elph` function extracted from
    YamboElectronPhononDB. You can make your own plotting functions.
    """
    fig, ax = plt.subplots(1,1)

    # Draw 2D BZ hexagonal borders
    ax.add_patch(BZ_hexagon(rlat))

    if plt_cbar:
        if 'cmap' in kwargs.keys(): color_map = plt.get_cmap(kwargs['cmap'])
        else:                       color_map = plt.get_cmap('viridis')

    # Reproduce plot also in adjacent BZs
    BZs = shifted_grids_2D(pts,rlat)
    for pts_s in BZs: plot=ax.scatter(pts_s[:,0],pts_s[:,1],c=data,**kwargs)

    if plt_cbar: cbar = fig.colorbar(plot)

    plt.gca().set_aspect('equal')

    lim = np.linalg.norm(rlat[1])
    lim = 0.85*np.linalg.norm(rlat[0])
    ax.set_xticks([])
    ax.set_yticks([])
    ticks = [-0.15,-0.1,-0.05,0,0.05,0.1,0.15]
    ticklabels = [ "{:.1f}".format(tick*2.*np.pi/bohr2ang) for tick in ticks ]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(ticklabels)
    ax.set_yticklabels(ticklabels)
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_xlabel(r'$q_x$ ($\mathring{A}^{-1}$)',fontsize=12)
    ax.set_ylabel(r'$q_y$ ($\mathring{A}^{-1}$)',fontsize=12)
    #cbar.ax.set_ylabel(cbarlabel,fontsize=12)

if __name__ == "__main__":
    """
    Main part of the script
    """
    # Which database to read
    elphdb        = 'lelphc'
    elphdb        = 'yambo'
    
    # In which momentum space to plot
    kspace_Plot    = True # k-space at fixed q 
    qspace_Plot    = True # q-space at fixed k

    # MoS2 details
    # * Mode A_1^\prime --> Index 7 (penultimate mode)
    # * v bands --> Index 0 [24], 1 [25]
    # * c bands --> Index 2 [26], 3 [27]

    # Customly chosen matrix elements
    i_v, i_vp = [0,1] #We choose the two valence bands (spin-orbit split at K/K')
    i_nu = 7          # A1^\prime phonon mode 
    iK_qspace = 34    # indx of K-pt in 12x12x1 BZ for q-space (phonon momenta)
    iK_kspace = 34    # indx of K-pt in 12x12x1 BZ for k-space (electron momenta)

    # Create "lattice" object from ns.db1 db
    ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')
    # Restrict plot two q=0 plane (no effect in true 2D)
    indx2D=(ylat.car_kpoints[:,2]==0.).nonzero()[0]
    # In yambo, k and q grids coincide
    ptgrid = ylat.car_kpoints[indx2D]

    # Create "elphon" object by reading:
    # (i) either the ndb.elph_gkkp* databases inside the yambo SAVE
    # (ii) or the ndb.elph database from lelphc
    if elphdb=='lelphc':
        
        elph = LetzElphElectronPhononDB(save_path+'/ndb.elph')
        print(elph)

        # Print info on how to use this class
        print(elph.__doc__)

        el_ph_mat_squared = elph.gkkp_sq[:,:,:,0,:,:]
        units = 'Ry'

    elif elphdb=='yambo':
        
        elph = YamboElectronPhononDB(ylat,folder_gkkp=save_path+'/SAVE',save=save_path+'/SAVE')
        print(elph)

        # Print info on how to use this class
        print(elph.__doc__)
    
        el_ph_mat_squared = elph.gkkp_sq
        units= 'Ha'
    
    else:
        print(f"elphdb string must be 'yambo' or 'lelphc', you have {elphdb}")
        exit()

    # We plot the average of \sqrt{|g_{vv'm}(k,q)|^2} over v and v', output in meV    
    # We select a specific Q-point to plot |g(K)| in kspace
    g_of_k = get_plottable( el_ph_mat_squared[iK_qspace,:,i_nu,:,:],i_v,i_vp,units=units )

    # We select a specific K-point to plot |g(Q)| in qspace
    g_of_q = get_plottable( el_ph_mat_squared[:,iK_kspace,i_nu,:,:],i_v,i_vp,units=units )
   
    # Plots are customisable as needed using matplotlib
    if kspace_Plot:

        plot_2D_elph(ptgrid,g_of_k,rlat=ylat.rlat,plt_cbar=True,\
                     marker='H',s=700,cmap='magma')
        plt.title(r'$|g(k)|$ (meV)',fontsize=8)
        plt.savefig('g_kspace.png',dpi=200)
        plt.show()

    if qspace_Plot:

        plot_2D_elph(ptgrid,g_of_q,rlat=ylat.rlat,plt_cbar=True,\
                     marker='H',s=700,cmap='magma')
        plt.title(r'$|g(q)|$ (meV)',fontsize=8)
        plt.savefig('g_qspace.png',dpi=200)
        plt.show()

