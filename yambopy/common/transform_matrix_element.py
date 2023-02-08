# This file is part of yambopy
# Author: FP
from yambopy import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from yambopy.units import ha2ev, ev2cm1, I
from yambopy.plot.plotting import add_fig_kwargs,BZ_hexagon,shifted_grids_2D

@add_fig_kwargs
def plot_BZ_2D(self,data,bzgrid_cc=None,plt_show=False,plt_cbar=False,**kwargs):
    """
    2D scatterplot in the BZ of the quantity A_{k}(i1,i2,i3,...).

    Any real quantity which is a function of only the k-grid may be supplied.
    The additional indices are user-specified.

    Example:
    - Exciton coefficients A_{k}(iq,ilambda,iv,ic) or A_{q}(ik,ilambda,iv,ic)
    - Electron-phonon matrix elements g_{k}(iq,imu,ib1,ib2) or g_{q}(ik,imu,ib1,ib2)
    - In general any k- or q-dependent quantity

    Inputs:
    - latticeDB object (lattice)
    - data do visualise    
    - BZ-grid in Cartesian coordinates (if not provided, uses lat.car_kpoints)

    Options:
    - if plt_show plot is shown
    - if plt_cbar colorbar is shown
    - kwargs example: marker='H', s=300, cmap='viridis', etc.

    Returns:
    - fig and ax objects

    NB: so far works for hexagonal systems. Can be improved to plot BZ planes at constant k_z for 3D systems and non-hexagonal cells.
    """
    rlat = lattice.rlat
    if bzgrid_cc is not None: 
        kpts = bz_grid_cc
    else:
        print("Using yambopy kpoints for plot") 
        kpts = self.car_kpoints

    # Input check
    if len(data)!=len(kpts):
        raise ValueError('Something wrong in data dimensions (%d data vs %d kpts)'%(len(data),len(kpts)))

    # Global plot stuff
    fig, ax = plt.subplots(1, 1)
    ax.add_patch(BZ_hexagon(rlat))

    if plt_cbar:
        if 'cmap' in kwargs.keys(): color_map = plt.get_cmap(kwargs['cmap'])
        else:                       color_map = plt.get_cmap('viridis')
    lim = 1.05*np.linalg.norm(rlat[0])
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)

    # Reproduce plot also in adjacent BZs
    BZs = shifted_grids_2D(kpts,rlat)
    for kpts_s in BZs: plot=ax.scatter(kpts_s[:,0],kpts_s[:,1],c=data,**kwargs)

    if plt_cbar: fig.colorbar(plot)

    plt.gca().set_aspect('equal')

    if plt_show: plt.show()
    else: print("Plot ready.\nYou can customise adding savefig, title, labels, text, show, etc...")

    return fig, ax

class ExpandMatrixElement():
    """
    This is a generic symmetry expansion script.

    Input: 
      - System symmetries from YamboLatticeDB
      - Array of matrix elements on k/q-grids to be expanded from IBZ to BZ
      - Instruction to expand over k or over q

    Output:
      - Array of expanded matrix elements in the specified grid

    Example:

      :: O_nmR{k}^S{q} =  <n R{k}|O|mR{k}-S{q}>

      :: Expansion over q is

      :: O_nmR{k}^S{q} =  O_nm S^-1R{k}^q     (no TR)
      :: O_nmR{k}^S{q} = [O_nm S^-1R{k}^q]^*  (TR)

      :: Expansion over k is

      :: O_nmR{k}^S{q} =  O_nm k^R^-1S{q}     (no TR)
      :: O_nmR{k}^S{q} = [O_nm k^R^-1S{q}]^*  (TR)

    """
    #def __init__(mats_ibz,syms,space='q',TR=False):
        
