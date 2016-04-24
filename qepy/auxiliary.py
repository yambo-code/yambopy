# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
import numpy as np
#try:
#    from matplotlib import pyplot as plt

def car_red(car,lat):
    """
    Convert cartesian coordinates to reduced
    """
    return np.array(map( lambda coord: np.linalg.solve(np.array(lat).T,coord), car))

def generate_path(klist,kinterval):
    """
    Generation of a path in reciprocal space by specifying a list of k-points
    Output in the format of quantum espresso == [ [kx, ky, kz, 1], ... ]
    Probably Henrique is going to suffer a hearth attack...
    """
    kout  = np.zeros([sum(kinterval)+1,4])
    kout[:,3] = 1
    #klist = np.array(klist)
    klabel = []
    kpoint = []
    for kline in klist:
      klabel.append(kline[1])
      kpoint.append(kline[0]) 
    kpoint = np.array(kpoint)
    io = 0
    for ik,interval in enumerate(kinterval):
      for ip in range(interval):
        kout[io,:3] = kpoint[ik] + float(ip)/interval*(kpoint[ik+1] - kpoint[ik])
        io = io + 1
    kout[io,:3] = kpoint[ik] + float(ip+1)/interval*(kpoint[ik+1] - kpoint[ik])
    return kout

def plot_qe_bands(prefix,folder,path,kint,nbands,nkpoints):
    """
    Plot of the band structure.
    """
    os.system('cd %s ; band_reading_qe %s' % (folder,prefix))
    bands = np.loadtxt('%s/%s.dat' % (folder,prefix))
    klabel, kpos = [], []
    fig = plt.subplots(figsize=(11,10))
    kpos.append(0)
    for line in path:
      klabel.append(line[1]) 
    for kl in kint:
      kpos.append(kpos[-1]+kl)
    plt.xticks( kpos, klabel)
    plt.ylabel('E (eV)')
    for ib in range(nbands):
      plt.plot(bands[ib*nkpoints:(ib+1)*nkpoints,0], bands[ib*nkpoints:(ib+1)*nkpoints,1], 'r-', lw=2)
    plt.show()
