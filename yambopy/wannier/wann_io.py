# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
# Author: Riccardo Reho 2023
# Adapted from wannier-berri

import numpy as np
import multiprocessing
import gc
from concurrent.futures import ProcessPoolExecutor
import tbmodels
from itertools import islice
from time import time
import fortio, scipy.io

# lambda function needed for reading Files
readstr = lambda F: "".join(c.decode('ascii') for c in F.read_record('c')).strip()

class W90_data():

    @property
    def n_neighb(self):
        return 0

    @property
    def NK(self):
        return self.data.shape[0]

    @property
    def NB(self):
        return self.data.shape[1 + self.n_neighb]

    @property
    def NNB(self):
        if self.n_neighb > 0:
            return self.data.shape[1]
        else:
            return 0

def convert(A):
    return np.array([l.split() for l in A], dtype=np.float64)

class MMN(W90_data):
    """
    MMN.data[ik, ib, m, n] = <u_{m,k}|u_{n,k+b}>
    """

    @property
    def n_neighb(self):
        return 1

    def __init__(self, seedname, npar=multiprocessing.cpu_count()):
        t0 = time()
        f_mmn_in = open(seedname + ".mmn", "r")
        f_mmn_in.readline()
        NB, NK, NNB = np.array(f_mmn_in.readline().split(), dtype=int)
        self.data = np.zeros((NK, NNB, NB, NB), dtype=np.complex128)
        block = 1 + self.NB * self.NB
        data = []
        headstring = []
        mult = 4
        # FIXME: npar = 0 does not work
        if npar > 0:
            pool = multiprocessing.Pool(npar)
        for j in range(0, NNB * NK, npar * mult):
            x = list(islice(f_mmn_in, int(block * npar * mult)))
            if len(x) == 0: break
            headstring += x[::block]
            y = [x[i * block + 1:(i + 1) * block] for i in range(npar * mult) if (i + 1) * block <= len(x)]
            if npar > 0:
                data += pool.map(convert, y)
            else:
                data += [convert(z) for z in y]

        if npar > 0:
            pool.close()
            pool.join()
        f_mmn_in.close()
        t1 = time()
        data = [d[:, 0] + 1j * d[:, 1] for d in data]
        self.data = np.array(data).reshape(self.NK, self.NNB, self.NB, self.NB).transpose((0, 1, 3, 2))
        headstring = np.array([s.split() for s in headstring], dtype=int).reshape(self.NK, self.NNB, 5)
        assert np.all(headstring[:, :, 0] - 1 == np.arange(self.NK)[:, None])
        self.neighbours = headstring[:, :, 1] - 1
        self.G = headstring[:, :, 2:]
        t2 = time()
        print("Time for MMN.__init__() : {} , read : {} , headstring {}".format(t2 - t0, t1 - t0, t2 - t1))

class AMN(W90_data):
    """
    The file seedname.amn contains the projection Amn(k).
    First line: a user comment, e.g., the date and time
    Second line: 3 integers: num_bands, num_kpts, num_wann
    Subsequently num_bandsxnum_wannxnum_kpts lines: 3 integers and 2 real numbers on each line. 
    The first two integers are the band index m and the projection index n, respectively. 
    The third integer specifies the k by giving the ordinal corresponding to its position in the list of k-points in seedname.win. 
    The real numbers are the real and imaginary parts, respectively, of the actual Amn(k).
    """
    def __init__(self, seedname, npar = multiprocessing.cpu_count()):
        t0 = time()
        f_amn_in = open(seedname + '.amn', "r")
        f_amn_in.readline()
        NB, NK, NW = np.array(f_amn_in.readline().split(),dtype=int) # number of bands, number of k-points, number of wannier
        f_amn_in.close()
        A_kmn = np.zeros((NB, NW, NK), dtype=np.complex128)
        print(NB)

        data = np.loadtxt(seedname + '.amn', dtype=np.float64, skiprows=2)
        A_kmn = data[:,3] + 1j*data[:,4]
        self.A_kmn = A_kmn.reshape(NK, NW, NB).transpose(0, 2, 1)

        t1 = time()
        print("Time for AMN.__init__() : {}".format(t1 - t0))
    
class EIG(W90_data):

    def __init__(self, seedname):
        data = np.loadtxt(seedname + ".eig")
        NB = int(round(data[:, 0].max()))
        NK = int(round(data[:, 1].max()))
        data = data.reshape(NK, NB, 3)
        assert np.linalg.norm(data[:, :, 0] - 1 - np.arange(NB)[None, :]) < 1e-15
        assert np.linalg.norm(data[:, :, 1] - 1 - np.arange(NK)[:, None]) < 1e-15
        self.data = data[:, :, 2]

class UXU(W90_data):
    """
    Read and setup uHu or uIu object.
    pw2wannier90 writes data_pw2w90[n, m, ib1, ib2, ik] = <u_{m,k+b1}|X|u_{n,k+b2}>
    in column-major order. (X = H for UHU, X = I for UIU.)
    Here, we read to have data[ik, ib1, ib2, m, n] = <u_{m,k+b1}|X|u_{n,k+b2}>.
    """

    @property
    def n_neighb(self):
        return 2

    def __init__(self, seedname='wannier90', formatted=False, suffix='uHu'):
        print("----------\n  {0}   \n---------".format(suffix))
        print('formatted == {}'.format(formatted))
        if formatted:
            f_uXu_in = open(seedname + "." + suffix, 'r')
            header = f_uXu_in.readline().strip()
            NB, NK, NNB = (int(x) for x in f_uXu_in.readline().split())
        else:
            f_uXu_in = FortranFileR(seedname + "." + suffix)
            header = readstr(f_uXu_in)
            NB, NK, NNB = f_uXu_in.read_record('i4')

        print("reading {}.{} : <{}>".format(seedname, suffix, header))

        self.data = np.zeros((NK, NNB, NNB, NB, NB), dtype=np.complex128)
        if formatted:
            tmp = np.array([f_uXu_in.readline().split() for i in range(NK * NNB * NNB * NB * NB)], dtype=np.float64)
            tmp_cplx = tmp[:, 0] + 1.j * tmp[:, 1]
            self.data = tmp_cplx.reshape(NK, NNB, NNB, NB, NB).transpose(0, 2, 1, 3, 4)
        else:
            for ik in range(NK):
                for ib2 in range(NNB):
                    for ib1 in range(NNB):
                        tmp = f_uXu_in.read_record('f8').reshape((2, NB, NB), order='F').transpose(2, 1, 0)
                        self.data[ik, ib1, ib2] = tmp[:, :, 0] + 1j * tmp[:, :, 1]
        print("----------\n {0} OK  \n---------\n".format(suffix))
        f_uXu_in.close()


class UHU(UXU):
    """
    UHU.data[ik, ib1, ib2, m, n] = <u_{m,k+b1}|H(k)|u_{n,k+b2}>
    """

    def __init__(self, seedname='wannier90', formatted=False):
        super().__init__(seedname=seedname, formatted=formatted, suffix='uHu')


class UIU(UXU):
    """
    UIU.data[ik, ib1, ib2, m, n] = <u_{m,k+b1}|u_{n,k+b2}>
    """

    def __init__(self, seedname='wannier90', formatted=False):
        super().__init__(seedname=seedname, formatted=formatted, suffix='uIu')

class FortranFileR(fortio.FortranFile):

    def __init__(self, filename):
        print("using fortio to read")
        try:
            super().__init__(filename, mode='r', header_dtype='uint32', auto_endian=True, check_file=True)
        except ValueError:
            print("File '{}' contains subrecords - using header_dtype='int32'".format(filename))
            super().__init__(filename, mode='r', header_dtype='int32', auto_endian=True, check_file=True)

class HR(W90_data):
    """
    HR.data[nrpts, num_wann,num_wann] = h_{Rmn}
    he second line states the number of Wannier functions num_wann. 
    The third line gives the number of Wigner-Seitz grid-points nrpts. 
    The next block of nrpts integers gives the degeneracy of each Wigner-Seitz grid point, with 15 entries per line. 
    Finally, the remaining num_wann2x nrpts lines each contain, respectively, 
    the components of the vector R in terms of the lattice vectors {Ai}, 
    the indices m and n, and the real and imaginary parts of the Hamiltonian matrix element Hmn(R) in the WF basis, e.g.,
    """

    def __init__(self, seedname, npar=multiprocessing.cpu_count()):
        t0 = time()
        file_path = seedname + "_hr.dat"
        ws_deg = []
        with open(file_path, "r") as f:
            f.readline()
            self.num_wann = int(f.readline())
            self.nrpts = int(f.readline())
            self.skiplines = int(np.ceil(self.nrpts/15))
            for i in range(0,int(np.ceil(self.nrpts/15))):
                ws_deg = np.append(ws_deg, f.readline().split())

        self.ws_deg = list(map(int,ws_deg))
        self.data = np.zeros((self.nrpts, self.num_wann,self.num_wann), dtype=np.complex128)
        data = np.loadtxt(file_path, dtype=np.float64, skiprows=self.skiplines+3)

        t1 = time()
        HR_mn = data[:,5] + 1j*data[:,6]
        iHR_mn = data[:,0:5]
        self.HR_mn = np.array(HR_mn).reshape(self.nrpts, self.num_wann, self.num_wann)
        self.iHR_mn = np.array(iHR_mn).reshape(self.nrpts,self.num_wann,self.num_wann,5)
        newhop = self.iHR_mn[:,:,:,0:3].reshape(self.nrpts*self.num_wann**2,3)
        self.hop = newhop[::self.num_wann**2]
        t2 = time()
        print("Time for HR.__init__() : {} , read : {} , headstring {}".format(t2 - t0, t1 - t0, t2 - t1))

class RMN(W90_data):
    """
    rmn.data[nrpts, num_wann,num_wann] = r_{mn}
    """

    def __init__(self, seedname, npar=multiprocessing.cpu_count()):
        t0 = time()
        f_hr_in = open(seedname + "_r.dat", "r")
        f_hr_in.readline()
        self.num_wann = int(f_hr_in.readline())
        self.nrpts = int(f_hr_in.readline())

        self.data = np.zeros((self.nrpts, self.num_wann,self.num_wann), dtype=np.complex128)
        block = self.num_wann**2
        data = []
        mult = 1
        # FIXME: npar = 0 does not work
        if npar > 0:
            pool = multiprocessing.Pool(npar)
        for j in range(0, self.nrpts, npar * mult):
            x = list(islice(f_hr_in, int(block * npar * mult)))
            if len(x) == 0: break
            y = [x[i * block :(i + 1) * block] for i in range(npar * mult) if (i + 1) * block <= len(x)]
            if npar > 0:
                data += pool.map(convert, y)
            else:
                data += [convert(z) for z in y]

        if npar > 0:
            pool.close()
            pool.join()
        f_hr_in.close()
        t1 = time()
        r_mn = [d[:, 5:11:2] + 1j * d[:, 6:11:2] for d in data]
        #print(r_mn[0])
        ir_mn = [d[:,0:5] for d in data]
        self.r_mn = np.array(r_mn).reshape(self.nrpts, self.num_wann, self.num_wann,3)
        self.ir_mn = np.array(ir_mn).reshape(self.nrpts,self.num_wann,self.num_wann,5)
        t2 = time()
        print("Time for RMN.__init__() : {} , read : {} , headstring {}".format(t2 - t0, t1 - t0, t2 - t1))

class NNKP():
    """
    read nnkp data file. 
    Store reciprocal lattice vectors, neighbouring points and b-vectors
    """
    def __init__(self, seedname):
        t0 = time()
        f_hr_in = open(seedname + ".nnkp", "r")
        real_lattice = []
        reciprocal_lattice = []
        k = []
        nkpoints = 0 
        nnkpts = 0
        b_grid = []
        data = []

        lines = f_hr_in.readlines()
        current_block = None

        for line in lines:
            if line.startswith('begin'):
                current_block = line.split()[1]
            elif line.startswith('end'):
                current_block = None
            else:
                if current_block == 'real_lattice':
                    real_lattice.append([np.float64(x) for x in line.split()])
                elif current_block == 'recip_lattice':
                    reciprocal_lattice.append([np.float64(x) for x in line.split()])
                elif current_block == 'kpoints':
                    if nkpoints == 0:
                        nkpoints = int(line.strip())
                    else:
                        k.append([np.float64(x) for x in line.split()])
                elif current_block == 'nnkpts':
                    if nnkpts == 0:
                        nnkpts = int(line.strip())
                    else:
                        data.append([np.float64(x) for x in line.split()])

        self.real_lattice = np.array(real_lattice)
        self.reciprocal_lattice = np.array(reciprocal_lattice)
        self.k = np.array(k)
        self.nkpoints = nkpoints 
        self.data = np.array(data)
        self.ik = self.data[:,0].astype(int)-1
        self.ikpb = self.data[:,1].astype(int)-1
        self.iG = self.data[:,2:5].astype(int)
        self.nnkpts = nnkpts

        Gvec = np.zeros((self.data.shape[0],3), dtype = np.longdouble)
        
        for ikkp in range(0, len(self.ik)):
            #Gvec[ikkp] = np.dot(reciprocal_lattice, self.iG[ikkp] )
            tmpb = self.k[self.ikpb[ikkp]] - self.k[self.ik[ikkp]] - self.iG[ikkp]
            b_grid.append(tmpb)

        self.b_grid = np.array(b_grid)

        t2 = time()
        print("Time for NNKP.__init__() : {}".format(t2 - t0))
