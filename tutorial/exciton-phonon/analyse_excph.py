"""
Simple script to plot and analyse 
exciton-phonon matrix elements calculated with yambo_ph
"""
import numpy as np
import matplotlib.pyplot as plt
from yambopy import YamboLatticeDB,YamboExcitonPhononDB
from yambopy.units import ha2ev,bohr2ang

def plot2D(excph,exc_in,exc_out,ph_in,plt_cbar=False,**kwargs):
    """ 
    Input indices are lists (even of one element) or 'all'
    
    Example: we can sum |G|^2 over all outgoing states and phonon modes 
             for the first bright exc_in state
    
    np.sum( X.excph_sq, axis=3 ) # sum over exc_out
    np.sum( data, axis=1 )    # sum over phonon modes
    np.sqrt(data[:,2]+data[:,3]) # sum over degenerate states (MoS2 A-bright exc_in)

    We construct F_q to plot
    """

    G_squared = excph.excph_sq
    G2plt = np.zeros(len(G_squared))

    if exc_in  == 'all': exc_in  = range(G_squared.shape[2])
    if exc_out == 'all': exc_out = range(G_squared.shape[3])
    if ph_in   == 'all': ph_in   = range(G_squared.shape[1])

    G_squared = G_squared[:, ph_in, :, :].sum(axis=(1))
    G_squared = G_squared[:, exc_in, :].sum(axis=(1))
    G_squared = G_squared[:, exc_out].sum(axis=(1))

    F_q = np.sqrt( G_squared )*ha2ev # Switch from Ha to eV

    # Do plot
    excph.plot_excph(F_q,plt_cbar=plt_cbar,**kwargs)
    lim = np.linalg.norm(excph.rlat[1])
    lim = 0.85*np.linalg.norm(excph.rlat[0])
    excph.ax.set_xticks([])
    excph.ax.set_yticks([])
    ticks = [-0.15,-0.1,-0.05,0,0.05,0.1,0.15]
    ticklabels = [ "{:.1f}".format(tick*2.*np.pi/bohr2ang) for tick in ticks ]
    excph.ax.set_xticks(ticks)
    excph.ax.set_yticks(ticks)
    excph.ax.set_xticklabels(ticklabels)
    excph.ax.set_yticklabels(ticklabels)
    excph.ax.set_xlim(-lim,lim)
    excph.ax.set_ylim(-lim,lim)
    excph.ax.set_xlabel('$q_x$ ($\mathring{A}^{-1}$)',fontsize=12)
    excph.ax.set_ylabel('$q_y$ ($\mathring{A}^{-1}$)',fontsize=12)
    excph.cbar.ax.set_ylabel('$|F(q)|$ eV',fontsize=12)

#
# User input
#

# Exciton in states
exc_in  = [2,3]     # A: 2,3 -- B: 6,7
exc_out = [0,1,2,3] # first 4 states (dispersion of triplet state and A)
ph_in  = 'all' 

# Paths of databases
path='.'
ns_db1 =f'{path}/SAVE/ns.db1'
ndb_exc=f'{path}/excph'

# Plot details
plt_fl='A_exc_2Dmap.png'
plt_title='A exciton'
size = 700

#
# Load data
#

# Read lattice and k-space info
ylat = YamboLatticeDB.from_db_file(filename=ns_db1)#,Expand=False)
print(ylat)

# Read exc-ph databases
X = YamboExcitonPhononDB(ylat,save_excph=ndb_exc)
print(X)

#
# Plot
#
plot2D(X,exc_in,exc_out,ph_in,plt_cbar=True,\
       marker='H',s=size,cmap='magma')
plt.title(plt_title,fontsize=8)
plt.savefig(plt_fl,dpi=200)
plt.show()
