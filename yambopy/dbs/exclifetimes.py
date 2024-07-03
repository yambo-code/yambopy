#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC, AMS, FP, RR
#
# This file is part of the yambopy project
#
from yambopy.units import *
from yambopy.plot.plotting import add_fig_kwargs,BZ_Wigner_Seitz
from yambopy.plot.bandstructure import *
from yambopy.lattice import replicate_red_kmesh, calculate_distances, get_path, car_red
from yambopy.tools.funcs import gaussian, lorentzian
from yambopy.tools.skw import SkwInterpolator
from yambopy.dbs.savedb import *
from yambopy.dbs.latticedb import *
from yambopy.dbs.electronsdb import *
from yambopy.dbs.qpdb import *
from yambopy.tools.skw import SkwInterpolator
from yambopy.dbs.dipolesdb import *

class ExcLifetimes():
        def __init__(self, yexcdb):
            self.yexcdb = yexcdb  # Initialize parent class

        def get_exciton_lifetimes(self, statelist = None, degen_step = 0.001, gauge = 'length', verbosity = 'False', pl_res = False):
            """get exciton lifetimes in 2D assuming the plane is on x,y
                TODO: T different than 0
            """
            if (statelist is None ):
                statelist = np.arange(1,self.yexcdb.nexcitons+1)
            Omega = self.yexcdb.lattice.lat_vol
            a1 = self.yexcdb.lattice.lat[0]
            a2 = self.yexcdb.lattice.lat[1]
            A = np.linalg.norm(np.cross(a1,a2))
            muS2 = 0 
            tmpE = np.ma.asarray(self.yexcdb.eigenvalues.real)
            tmpI = np.ma.asarray(self.yexcdb.l_residual * self.yexcdb.r_residual)
            excE = sorted(tmpE)
            excI = tmpI

            tau0_tot = np.zeros(len(statelist))
            merged_states = np.empty(len(statelist), dtype='object')
            # q0 norm factor 
            q0_norm = 1e-5

            for l, state in enumerate(statelist):
                state_internal = state 

                # find states within degen_step window from state
                mesk = np.logical_and(excE>=excE[state_internal]-degen_step,excE<=excE[state_internal]+degen_step)
                states_idx = np.where(mesk==True)[0]
                #compute gamma
                gamma0 = 0
                for i, st in enumerate(states_idx):
                    ES = excE[st]/ha2ev
                    if (gauge=='length'):
                        muS2 = excI[st]/(q0_norm**2) # if you inspect the Yambo code you might expect another 1/((2*np.pi)**3) but I think that d3k_factor/((2np.pi)**3) is actually 1/Omega
                    elif (gauge == 'velocity'):
                        muS2 = excI[st]/(ES**2)
                    gg = 4.*np.pi*ES*(muS2/self.yexcdb.nkpoints)/(A*speed_of_light)
                    gamma0 += gg
                #compute tau
                tau0_tot[l] = autime2s/(gamma0.real)
                merged_states[l] = '{}<->{}'.format(min(states_idx)+1,max(states_idx)+1)

            if (verbosity):
                for l,state in enumerate(statelist):
                    print(' Exciton radiative lifetime')
                    print('#')
                    #print(f'# Effective mass = {m_e}\n')
                    print(f"# Global Gauge = {gauge}")
                    print(f'Energy degeneration step = {degen_step} eV')
                    print('#')
                    print(f"# {'State':>5} {'Energy (eV)':>15} {'$tau_0$ (s)':>20} {'Merged states':>20}")
                    print(f"{state:>5} {excE[state]:>15} {tau0_tot[l]:>20} {merged_states[l]:>20}")
            return excE, tau0_tot, merged_states
