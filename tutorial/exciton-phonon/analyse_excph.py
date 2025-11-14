"""
Simple script to plot and analyse
exciton-phonon matrix elements calculated with yambopy
"""

import numpy as np
import matplotlib.pyplot as plt
from yambopy import YamboLatticeDB
from yambopy.plot.plotting import shifted_grids_2D,BZ_hexagon
from yambopy.units import ha2ev,bohr2ang

def GET_G2_to_plot(G_squared, exc_in, exc_out, ph_in):
    if exc_in  == 'all': exc_in  = range(G_squared.shape[2])
    if exc_out == 'all': exc_out = range(G_squared.shape[3])
    if ph_in   == 'all': ph_in   = range(G_squared.shape[1])

    G_squared = G_squared[:, ph_in, :, :].sum(axis=(1))
    G_squared = G_squared[:, exc_in, :].sum(axis=(1))
    G_squared = G_squared[:, exc_out].sum(axis=(1))

    F_q = np.sqrt( G_squared )*ha2ev
    return F_q
    
def plot_2D_excph(qpts,data,rlat=None,plt_cbar=False,**kwargs):

    fig, ax = plt.subplots(1,1)
    
    # Draw 2D BZ hexagonal borders
    ax.add_patch(BZ_hexagon(rlat))

    if plt_cbar:
        if 'cmap' in kwargs.keys(): color_map = plt.get_cmap(kwargs['cmap'])
        else:                       color_map = plt.get_cmap('viridis')

    # Reproduce plot also in adjacent BZs
    BZs = shifted_grids_2D(qpts,rlat)
    for qpts_s in BZs: plot=ax.scatter(qpts_s[:,0],qpts_s[:,1],c=data,**kwargs)

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
    ax.set_xlabel('$q_x$ ($\mathring{A}^{-1}$)',fontsize=12)
    ax.set_ylabel('$q_y$ ($\mathring{A}^{-1}$)',fontsize=12)
    cbar.ax.set_ylabel('$|\mathcal{G}(q)|$ eV',fontsize=12)

#
# User input
#
path = '3D_hBN'
path = '1L_MoS2'

# Exciton in states
exc_in  = [2,3]  # First bright peak
exc_out = [0,1,2,3]
ph_in  = 'all' 
# Paths of databases
ns_db1 =f'{path}/SAVE/ns.db1'
ns_ypy = 'MoS2_Ex-ph.npy'
# output file
plt_fl="MoS2_A_exc_2Dmag.png"

plt_title=f'Bright exciton | {ph_in} phonons, {exc_out} excitons'

#
# Load data
#

# Read lattice and k-space info
ylat = YamboLatticeDB.from_db_file(filename=ns_db1)
print(ylat)

# Restrict plot two q=0 plane (no effect in true 2D)
indx2D=(ylat.car_kpoints[:,2]==0.).nonzero()[0]

# Load exc-ph database
X_py = np.load(ns_ypy)
#X_py = X_py[:,:,:4,:12]
G_squared = np.abs(X_py)**2.

# Prepare quantity to plot
G2_to_plot = GET_G2_to_plot(G_squared,exc_in,exc_out,ph_in)[indx2D]
qgrid = ylat.car_kpoints[indx2D]

# Plot
plot_2D_excph(qgrid,G2_to_plot,rlat=ylat.rlat,plt_cbar=True,\
              marker='H',s=700,cmap='magma')
plt.title(plt_title,fontsize=8)
plt.savefig(plt_fl,dpi=200)
plt.show()
