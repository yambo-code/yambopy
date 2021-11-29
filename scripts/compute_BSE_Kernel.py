# Copyright (C) 2021 Davide Romanin
# All rights reserved.
#
# This file is part of yambopy
#
import os
import numpy as np
from yambopy import *


# Lattice information
lat  = YamboLatticeDB.from_db_file(filename='./ns.db1')

Q = range(1,52)

f = open('Singlet_exchange.dat', 'w')
f.write("#q Re[i] Im[i], with i running over exciton states\n")
f.close()
g = open('Singlet_direct.dat', 'w')
g.write("#q Re[i] Im[i], with i running over exciton states\n")
g.close()
h = open('Triplet_direct.dat', 'w')
h.write("#q Re[i] Im[i], with i running over exciton states\n")
h.close()

for q in Q:
    f = open('Singlet_exchange.dat', 'a')
    g = open('Singlet_direct.dat', 'a')
    h = open('Triplet_direct.dat', 'a')
    # Singlet Exciton database read from db file
    ySing = YamboExcitonDB.from_db_file(lat,filename='ndb.BS_diago_Q'+str(q),folder='./BSE_Singlet')
    # Triplet Exciton database read from db file
    yTrip = YamboExcitonDB.from_db_file(lat,filename='ndb.BS_diago_Q'+str(q),folder='./BSE_Triplet')

    # BSE Exchange Interaction read from db file
    yExchange = YamboBSEKernelDB.from_db_file(lat,filename='ndb.BS_PAR_Q'+str(q),folder='./BSE_Exchange_only')
    # BSE Direct   Interaction read from db file
    yDirect   = YamboBSEKernelDB.from_db_file(lat,filename='ndb.BS_PAR_Q'+str(q),folder='./BSE_Triplet')


    # List of excitons to be analised
    states = [1,2]

    Vs = YamboBSEKernelDB.get_kernel_value(q,states,yExchange,ySing)
    Ws = YamboBSEKernelDB.get_kernel_value(q,states,yDirect,ySing)
    np.savetxt(f, Vs)
    np.savetxt(g, Ws)


    Wt = YamboBSEKernelDB.get_kernel_value(q,states,yDirect,yTrip)
    np.savetxt(h, Wt)

    h.close()
    g.close()
    f.close()
    print("Q=%i : Done\n" %(q))
