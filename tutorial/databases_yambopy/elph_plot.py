"""
Tutorial for YamboElectronPhononDB.

Electron-phonon matrix element reading and plotting

EDIT the path below to point to the yambo SAVE folder.
"""
save_path='ELPH_saves/YAMBO_saves'
from yambopy import *
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    """
    Main part of the script
    """
    Kspace_Plot    = True
    Qspace_Plot    = False

    # Customly chosen matrix elements
    i_n, i_m = [3,4] #i_n =3 -> valence band, i_n = 4 -> conduction band in 2D-hBN
    i_nu = 3 # LA phonon mode in 2D-hBN at K (ZO mode at Gamma); LO, TO modes are i_nu=4,5
    i_q = 143 # This is the K-point in the hexagonal BZ
    i_k = 143

    #                    #
    # Start Yambopy part #
    #                    #

    # Create "lattice" object by reading the ns.db1 database inside the yambo SAVE
    ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')

    # Create "elphon" object by reading the ndb.elph_gkkp* databases inside the yambo SAVE
    yelph = YamboElectronPhononDB(ylat,folder_gkkp=save_path+'/SAVE',save=save_path+'/SAVE')
    print(yelph)

    # Print info on how to use this class
    print(yelph.__doc__)

    # We select a specific Q-point to plot |g(K)| in kspace
    g_of_k = np.abs(yelph.gkkp[i_q,:,i_nu,i_n,i_m])

    # We select a specific K-point to plot |g(Q)| in qspace
    g_of_q = np.abs(yelph.gkkp[:,i_k,i_nu,i_n,i_m])
    
    # Plots are customisable as needed using matplotlib
    if Kspace_Plot:
        yelph.plot_elph(g_of_k,s=100,plt_cbar=False,marker='H',cmap='viridis')
        yelph.ax.set_title('|g(k)| (Hartree)')
        yelph.ax.set_xlabel('k_x')
        yelph.ax.set_ylabel('k_y')
        plt.show()

    if Qspace_Plot:
        yelph.plot_elph(g_of_q,s=100,plt_cbar=False,marker='H',cmap='viridis')
        yelph.ax.set_title('|g(q)| (Hartree)')
        yelph.ax.set_xlabel('q_x')
        yelph.ax.set_ylabel('q_y')
        plt.show()

    #                   #
    # End Yambopy part. #
    #                   # 
