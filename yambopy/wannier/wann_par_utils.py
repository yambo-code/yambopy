import numpy as np
from yambopy.dbs.excitondb import *
from yambopy.dbs.bsekerneldb import *


def get_exc_overlap_ttp_par(argslist):
    l, lp, iq, ikq, ib, eigvec, h2peigvec_vck, BSE_table, kmpgrid, nv, bse_nv, dimbse = argslist
    Mssp_ttp = 0
    ik, iv, ic = BSE_table[l]
    ikp, ivp, icp = BSE_table[lp]
    iqpb = kmpgrid.qpb_grid_table[iq, ib][1]

    for it in range(dimbse):
        for itp in range(dimbse):
            ik = BSE_table[l][0]
            iv = BSE_table[l][1] 
            ic = BSE_table[l][2] 
            ikp = BSE_table[lp][0]
            ivp = BSE_table[lp][1] 
            icp = BSE_table[lp][2] 
            iqpb = kmpgrid.qpb_grid_table[iq, ib][1]
            ikmq = kmpgrid.kmq_grid_table[ik,iq][1]
            ikpbover2 = kmpgrid.kpbover2_grid_table[ik, ib][1]
            ikmqmbover2 = kmpgrid.kmqmbover2_grid_table[ik, iq, ib][1]
            Mssp_ttp += np.conjugate(h2peigvec_vck[ikq,l,bse_nv-nv+iv,ic-nv,ik])*h2peigvec_vck[iqpb,lp,bse_nv-nv+ivp,icp-nv, ikpbover2]*\
                            np.vdot(eigvec[ik,:, ic], eigvec[ikpbover2,:, icp])*np.vdot(eigvec[ikmqmbover2,:,ivp], eigvec[ikmq,:,iv]) 
    return (l, lp, iq, ib, Mssp_ttp)


def process_file(args):
    idx, exc_db_file, data_dict = args
    # Unpacking data necessary for processing
    latdb, kernel_path, kpoints_indexes, HA2EV, BSE_table, kplusq_table, kminusq_table_yambo, eigv, f_kn = data_dict.values()

    yexc_atk = YamboExcitonDB.from_db_file(latdb, filename=exc_db_file)
    kernel_db = YamboBSEKernelDB.from_db_file(latdb, folder=f'{kernel_path}', Qpt=kpoints_indexes[idx]+1)
    aux_t = np.lexsort((yexc_atk.table[:,2], yexc_atk.table[:,1],yexc_atk.table[:,0]))
    K_ttp = kernel_db.kernel[aux_t][:,aux_t]  
    H2P_local = np.zeros((len(BSE_table), len(BSE_table)), dtype=np.complex128)

    for t in range(len(BSE_table)):
        ik, iv, ic = BSE_table[t]
        for tp in range(len(BSE_table)):
            ikp, ivp, icp = BSE_table[tp]
            ikplusq = kplusq_table[ik, kpoints_indexes[idx],1]
            ikminusq = kminusq_table_yambo[ik, kpoints_indexes[idx]]
            ikpminusq = kminusq_table_yambo[ikp, kpoints_indexes[idx]]
            K = -(K_ttp[t, tp]) * HA2EV
            deltaE = eigv[ik, ic] - eigv[ikpminusq, iv] if (ik == ikp and icp == ic and ivp == iv) else 0.0
            occupation_diff = -f_kn[ikpminusq, ivp] + f_kn[ikp, icp]
            element_value = deltaE + occupation_diff * K
            H2P_local[t, tp] = element_value
    return idx, H2P_local

