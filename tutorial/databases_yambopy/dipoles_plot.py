"""
Tutorial for YamboElectronPhononDB.

Electron-phonon matrix element reading and plotting

EDIT the path below to point to the yambo SAVE folder.
"""
save_path='BSE_saves/YAMBO_saves'
dipoles_path='BSE_saves/BSE_databases'
from yambopy import *
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    Kspace_Plot     = True
    Absorption_Plot = False

    # Customly chosen  matrix elements
    i_c, i_v = [4,3] # 3-> top valence, 4-> bottom conduction
    i_x, i_y = [0,1] # Cartesian direction: 0 -> x, 1 -> y, 2 -> z

    #                    #
    # Start Yambopy part #
    #                    #

    # Create "lattice" object by reading the ns.db1 database inside the yambo SAVE
    ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')

    # Read dipole matrix elements
    ydip = YamboDipolesDB(ylat,save=dipoles_path,filename='ndb.dipoles')
    print(ydip)

    # Read electron energies
    yel = YamboElectronsDB(ylat,save=save_path+'/SAVE')
    print(yel)
    
    # Plot dipoles in k-space (modulus, summed over x,y,z)
    dip_of_k = np.sqrt( np.sum( np.abs(ydip.dipoles[:,:,i_c,i_v])**2., axis=1) )
    #dip_of_k = np.sqrt( np.abs(ydip.dipoles[:,i_x,i_c,i_v])**2.+np.abs(ydip.dipoles[:,i_y,i_c,i_v])**2. ) # sum over x,y plane only
    ydip.plot_dipoles(dip_of_k,s=100,plt_cbar=False,marker='H',cmap='viridis')

    # Get independent-particles absorption spectrum
    w, eps2 = ydip.ip_eps2(yel,ntot_dip=10,broad=0.01,emin=0.,emax=7.,esteps=1000)

    #                   #
    # End Yambopy part. #
    #                   #
    
    if Kspace_Plot:
        # Plot is customisable as needed using matplotlib
        ydip.ax.set_title('|d(k)| (a.u.)')
        ydip.ax.set_xlabel('k_x')
        ydip.ax.set_ylabel('k_y')
        plt.show()

    if Absorption_Plot:
        # Plot is customisable as needed using matplotlib
        plt.title('IP absorption spectrum')
        plt.xlim(0,7)
        plt.xlabel('Energy (eV)')
        plt.ylim(0,np.max(eps2)*1.1)
        plt.plot(w,eps2,color='blue')
        plt.show()

