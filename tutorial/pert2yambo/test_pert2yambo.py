import os 
from netCDF4 import Dataset
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import h5py
from itertools import product
import pickle
import yaml 
from yambopy import *
from yambopy.pert2yambo.RTDB import *

save_path="./Yambo/tmp/gaas.save/"
ndb = "./ndb.RT_carriers"
cdyna_e = "./cdyna-elec/gaas_cdyna.h5"
cdyna_h = "./cdyna-hole/gaas_cdyna.h5"
teth5_e = "./cdyna-elec/gaas_tet.h5"
teth5_h = "./cdyna-hole/gaas_tet.h5"
yamlfile_e = "./cdyna-elec/gaas_dynamics-run.yml"
yamlfile_h = "./cdyna-hole/gaas_dynamics-run.yml"


#Debug=False
Debug=True


def generate_grid(grid):
    dx=1.0/float(grid[0])
    dy=1.0/float(grid[1])
    dz=1.0/float(grid[2])
    nkpt=np.prod(grid)
    k_grid=zeros([nkpt,3],float)
    ic=0
    for kx in range(0, grid[0]):
        for ky in range(0, grid[1]):
            for kz in range(0, grid[2]):
                k_grid[ic,:]=[float(kx)*dx,float(ky)*dy,float(kz)*dz]
                ic=ic+1
    return k_grid


# Return the integer distance for periodic grids
# defined by kgrid
def periodic_dist(ikpt1,ikpt2,kgrid):
    idist=ikpt1-ikpt2
    # distances are between -kgrid/2, and kgrid/2
    for idx in range(3):
        if idist[idx]<-int(kgrid[idx]/2):
            idist[idx]=idist[idx]+kgrid[idx]
        if idist[idx]>=int(kgrid[idx]/2):
            idist[idx]=idist[idx]-kgrid[idx]
    return idist

dynoccups = dynamic_occupations(tmp_out="./tmp",dyn_yamlfile=yamlfile_e,cdynafile=cdyna_e,teth5file=teth5_e,ndbfile=save_path+'/SAVE/'+ndb)
dynoccups.pert_grid_reduced()

ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')
y_k_grid=ylat.k_grid
p_k_grid=dynoccups.kpts_grid

### Generate fake grids to check the code ##############
if Debug:
    y_grid=np.array([4,1,1])
    p_grid=np.array([12,1,1])
    
    y_k_grid=y_grid
    p_k_grid=p_grid

#######################################################


print("Perturbo k-grid ",p_k_grid)
print("Yambo k-grid ",y_k_grid)

if any(p_k_grid%y_k_grid != 0):
    print("Incompatible k-grid!! ")
    sys.exit(0)
else:
    print("Compatible k-grids found")
    
grid_ratio=np.int32(p_k_grid/y_k_grid)
print("Grid ratio ",grid_ratio)
if(any(grid_ratio%2 ==0)):
   print("Please use an odd grid ratio! ")
   sys.exit(0)

dynoccups.read_pert_num_kpts()
if Debug:
    n_kpt_pert=np.prod(p_k_grid)
else:
    n_kpt_pert=dynoccups.num_pert_kpts
print("Number of k-points in perturbo : ",n_kpt_pert)

if Debug:
    pert_kpts=generate_grid(p_k_grid)
else:
    pert_kpts=dynoccups.read_perturbo_kpts()


pert_kpts=make_kpositive(pert_kpts.tolist())
small_q=np.zeros(3,dtype=float)
small_q=1.0/p_k_grid

pert_ikpt=[np.int32(np.rint(kpt/small_q)) for kpt in pert_kpts]
if Debug:
    with open('perturbo_ik_bz.pts', 'w') as f:
        f.write("#Perturbo ik-points in the BZ\n")
        for ikpt in pert_ikpt:
            f.write(str(ikpt)+'\n')


yambo_ikpt=[np.int32(np.rint(kpt/small_q)) for kpt in ylat.red_kpoints]
if Debug:
    with open('yambo_ik_bz.pts', 'w') as f:
        f.write("#Yambo ik-points in the BZ\n")
        for ikpt in yambo_ikpt:
            f.write(str(ikpt)+'\n')

yambo_ibz_kpts=ylat.get_ibz_kpoints(units='red')
print("Number of IBZ k-points in Yabmo : ",len(yambo_ibz_kpts))

yambo_ikpt_ibz=[np.int32(np.rint(kpt/small_q)) for kpt in yambo_ibz_kpts]
if Debug:
    with open('yambo_ik_ibz.pts', 'w') as f:
        f.write("#Yambo ik-points in the IBZ\n")
        for ikpt in yambo_ikpt_ibz:
            f.write(str(ikpt)+'\n')

#search perturbo neighboars for each yambo point
yneighboars=[]
for iyk in tqdm(range(len(yambo_ikpt_ibz))):
    ykpt=yambo_ikpt_ibz[iyk]
    ylist=[]
    for ipk,pkpt in enumerate(pert_ikpt):
        dist=periodic_dist(ykpt,pkpt,p_k_grid)
        if all(abs(dist)<=grid_ratio/2):
            ylist.append(ipk)
    yneighboars.append(ylist)

#Average number of neighboars for each y-kpts
ave_n=0
max_n=-1
min_n=1000000
for ylist in yneighboars:
    y_n=len(ylist)
    ave_n=ave_n+y_n
    if max_n < y_n:
        max_n=y_n
    if min_n > y_n:
        min_n=y_n
ave_n=ave_n/len(yambo_ikpt)
print("Average number of neighboards : ",ave_n)
print("Max/Min number of neighboards : ",max_n,min_n)

dynoccups.get_vcb_indices()
# dynoccups.parse_bands_from_yaml()
# dynoccups.get_files()


# Copy bare occupation for all k-points

# Update occupation only for the points included in the dynamics
