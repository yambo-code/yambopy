"""
Tutorial for YamboLatticeDB.

We show a specific functionality, namely kpoint reading and expansion.

EDIT the path below to point to the yambo SAVE folder.
"""
save_path='BSE_saves/YAMBO_saves'
from yambopy import *
from yambopy.plot.plotting import BZ_Wigner_Seitz
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    """ 
    Main part of the script
    """
    Cartesian_Plot = True
    Crystal_Plot   = False
    Symmetry_Plot  = False
    i_k = 8 #Custom kpoint index in the IBZ. To be expanded in the full BZ.

    #                    #
    # Start Yambopy part #
    #                    #

    # Create "lattice" object by reading the ns.db1 database inside the yambo SAVE
    ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')
    print(ylat)

    # Read it again without expansion of the kpts from IBZ to BZ
    ylat_noexp = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1',Expand=False)

    # Apply symmetries: we select one point in the IBZ and manually expand it by applying
    #                   all symmetries in Cartesian coordinates
    kpt = ylat_noexp.car_kpoints[i_k]
    klist = np.array([ np.matmul(symmetry, kpt) for symmetry in ylat_noexp.sym_car  ])

    #                                                                #
    # End Yambopy part. After this is just plotting with matplotlib. #
    #                                                                #

    if Cartesian_Plot:
        # Do a 2D scatterplot of the kpoints in Cartesian coordinates
        fig = plt.figure(figsize=(9,9))
        ax = plt.gca()
        ax.set_title('Cartesian coordinates')

        ## Add BZ borders
        ax.add_patch(BZ_Wigner_Seitz(ylat,color='black',linewidth=1.))

        ## Plot with "nice" layout
        ax.set_aspect('equal')
        ax.scatter(ylat.car_kpoints[:,0],ylat.car_kpoints[:,1],marker='H',s=200,color='teal',\
                   linewidth=0.5,edgecolors='black',label='expanded')
        ax.scatter(ylat_noexp.car_kpoints[:,0],ylat_noexp.car_kpoints[:,1],marker='h',s=100,color='orange',\
                   linewidth=0.5,edgecolors='black',label='unexpanded')

        ## Explicitly show kpt indices
        for i_k,kpt in enumerate(ylat.car_kpoints): 
            kx,ky = kpt[0],kpt[1]
            ax.annotate(i_k, (kx,ky), color='teal', xytext=(kx+0.003,ky+0.005))
        for i_k,kpt in enumerate(ylat_noexp.car_kpoints): 
            kx,ky = kpt[0],kpt[1]
            ax.annotate(i_k, (kx,ky), color='orange', xytext=(kx-0.005,ky+0.005))

        plt.legend()
        plt.show()

    if Crystal_Plot:
        # Do a 2D scatterplot of the kpoints in crystal coordinates
        fig = plt.figure(figsize=(9,9))
        ax = plt.gca()
        ax.set_title('Crystal coordinates')

        ax.set_aspect('equal')
        ax.scatter(ylat.red_kpoints[:,0],ylat.red_kpoints[:,1],marker='H',s=200,color='teal',\
                   linewidth=0.5,edgecolors='black',label='expanded')
        ax.scatter(ylat_noexp.red_kpoints[:,0],ylat_noexp.red_kpoints[:,1],marker='h',s=100,color='orange',\
                   linewidth=0.5,edgecolors='black',label='unexpanded')

        plt.legend()
        plt.show()

    if Symmetry_Plot:
        # Plot star of kpt
        fig = plt.figure(figsize=(9,9))
        ax = plt.gca()
        ax.set_title('Kpoint expansion')
        ax.add_patch(BZ_hexagon(ylat.rlat,color='black',linewidth=1.))
        ax.set_aspect('equal')
        ax.scatter(ylat.car_kpoints[:,0],ylat.car_kpoints[:,1],marker='H',s=200,color='teal',\
                   linewidth=0.5,edgecolors='black',label='expanded')

        ax.scatter(klist[:,0],klist[:,1],marker='h',s=100,color='orange',\
                   linewidth=0.5,edgecolors='black',label='Custom kpoint')

        plt.legend()
        plt.show()
