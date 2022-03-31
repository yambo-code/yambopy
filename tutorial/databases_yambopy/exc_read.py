"""
Tutorial for YamboExcitonDB.

Reading excitonic quantities from ndb.BS_diago_Q1

EDIT the path below to point to the yambo SAVE folder.
"""
save_path='BSE_saves/YAMBO_saves'
bse_path ='BSE_saves/BSE_databases'
from yambopy import *
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    #Customly chosen Q-point
    iQ=0 # 0-> Gamma point, i.e., optical absorption limit

    #                    #
    # Start Yambopy part #
    #                    #

    # Create "lattice" object by reading the ns.db1 database inside the yambo SAVE
    ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')

    # Read exciton data at Q=iQ
    yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q1')
    print(yexc)

    # Eigenvalues (exciton energies) and intensities (residuals squared)
    print('\nEquivalent of ypp -e s -b 1: ')
    energies    = yexc.eigenvalues.real
    intensities = yexc.get_intensities().real 
    for i_exc in range(10): print(i_exc+1,' %2.4f'%energies[i_exc],' %2.4f'%intensities[i_exc])
    print('...\n ')

    # Eigenvectors (exciton wave functions)
    print('Eigenvector shape (number of excitons, number of transitions): ')
    print(yexc.eigenvectors.shape,'\n ') 

    # Table (from transition basis to single particle basis)
    print('Transition index t = (kvc) -> Single-particle indices k, v, c ')
    for it, t in enumerate(yexc.table[:10]): print(it, ' -> ', t[0],t[1],t[2])
    print('...')

    #                    #
    # Start Yambopy part #
    #                    #
