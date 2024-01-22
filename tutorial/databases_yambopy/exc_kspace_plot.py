"""
Tutorial for YamboExcitonDB.

Plotting exciton wavefunction components in the BZ

EDIT the path below to point to the yambo SAVE folder.
"""
save_path='BSE_saves/YAMBO_saves'
bse_path ='BSE_saves/BSE_databases'
from yambopy import *
from qepy import *
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    Kspace_Plot = False
    Bands_Plot  = False
    Bands_Plot_Interpolate = True
    
    # Customly chosen Q-point
    iQ=0 # 0-> Gamma point, i.e., optical absorption limit
    # States to be merged together (because they are degenerate)
    #
    # You may try the following states:
    # 
    # [1,2], [3,4], [5], [6,7]
    #
    states = [1,2]

    #                    #
    # Start Yambopy part #
    #                    #

    # Create "lattice" object by reading the ns.db1 database inside the yambo SAVE
    ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')

    # Read exciton data at Q=iQ
    yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q1')

    # Plot of exciton weights in k-space
    if Kspace_Plot:
        fig = plt.figure(figsize=(6,6))
        ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])
        yexc.plot_exciton_2D_ax(ax,states,mode='hexagon',limfactor=0.8,scale= 320)
        plt.show()

    # Plot on top of the band structure
    
    ## [1.] Define path in crystal coordinates using class Path

    npoints = 20
    path = Path([ [[  0.0,  0.0,  0.0],'$\Gamma$'],
                  [[  0.5,  0.0,  0.0],'M'],
                  [[1./3.,1./3.,  0.0],'K'],
                  [[  0.0,  0.0,  0.0],'$\Gamma$']], 
                  [int(npoints*2),int(npoints),int(sqrt(5)*npoints)] )

    ## [2.] Read electron energies
    yel = YamboSaveDB.from_db_file(folder=save_path+'/SAVE')
    

    ## [3.A] Plot without interpolating the values
    if Bands_Plot:
        fig = plt.figure(figsize=(4,6))
        ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

        exc_on_bands = yexc.get_exciton_bs(yel,path,states,size=1.0)
        exc_on_bands.plot_ax(ax,c_bands='grey',c_weights='red')

        ax.set_ylim(-7.5,12.)
        plt.show()

    ## [3.B] Interpolate the values
    if Bands_Plot_Interpolate:
        fig = plt.figure(figsize=(4,6))
        ax  = fig.add_axes( [ 0.15, 0.15, 0.80, 0.80 ])

        # In case of problems with the interpolation, try to increase lpratio
        exc_on_bands = yexc.interpolate(yel,path,states,lpratio=10,f=None,size=0.5,verbose=True)
        exc_on_bands.plot_ax(ax,c_bands='grey',c_weights='red',alpha_weights=0.5)

        ax.set_ylim(-7.5,12.)
        plt.show()
    
    #                  #
    # End Yambopy part #
    #                  #

