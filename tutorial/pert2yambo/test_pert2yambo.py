import os 
from netCDF4 import Dataset
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import h5py
from itertools import product
import shutil
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


Debug=False
# Debug=True

def copy_occupation(yneighboars,occups,occ_key,delta_f,otype):
    if otype!='e' and otype!='h':
        print("Error occupation can be only electron or hole ")
        sys.exit(0)

#Copy Electron/Hole Occupation
    n_bands=len(occups.pert_bands)
    my_occ=np.zeros(n_bands,float)
    for y_ikp,ylist in enumerate(yneighboars):
        my_occ[:]=0.0
        for k_idx in ylist:
            for idx,bnd in enumerate(occups.pert_bands):
                my_occ[idx]=my_occ[idx]+occups.occupation[occ_key][k_idx,idx]
                # I have to remove the bare occupation
                if otype=='h':
                    my_occ[idx]=my_occ[idx]-1.0
        
        #Copy occupation in a Yambo shape array delta_f
        for idx,bnd in enumerate(occups.pert_bands):
            delta_f[y_ikp,bnd]=delta_f[y_ikp,bnd]+my_occ[idx]

def generate_grid(grid):
    dx=1.0/float(grid[0])
    dy=1.0/float(grid[1])
    dz=1.0/float(grid[2])
    nkpt=np.prod(grid)
    k_grid=np.zeros([nkpt,3],float)
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
        if kgrid[idx]==1:
            idist[idx]=0
        else:
            if idist[idx]<-int(kgrid[idx]/2):
                idist[idx]=idist[idx]+kgrid[idx]
            if idist[idx]>=int(kgrid[idx]/2):
                idist[idx]=idist[idx]-kgrid[idx]
    return idist


def find_neighboars(y_ikpt,p_ikpt,p_kgrid,py_grid_ratio,Debug=False):
    y_neighboars=[]
    for iyk in tqdm(range(len(y_ikpt))):
        ykpt=y_ikpt[iyk]
        ylist=[]
        if Debug:
            i_nn=0
        for ipk,pkpt in enumerate(p_ikpt):
            dist=periodic_dist(ykpt,pkpt,p_kgrid)
            if all(abs(dist)<=py_grid_ratio/2):
                if Debug:
                    i_nn=i_nn+1
                ylist.append(ipk)
        if Debug:
            print("Number of ",iky,"  = ",i_nn)
        y_neighboars.append(ylist)

    # Report the average number of neighboars for each y-kpts
    ave_n=0
    max_n=-1
    min_n=1000000
    for ylist in y_neighboars:
        y_n=len(ylist)
        ave_n=ave_n+y_n
        if max_n < y_n:
            max_n=y_n
        if min_n > y_n:
            min_n=y_n
    ave_n=ave_n/len(yambo_ikpt_ibz)
    print("Average number of neighboards : ",ave_n)
    print("Max/Min number of neighboards : ",max_n,min_n)
    return y_neighboars



# read electrons
elec_occups = dynamic_occupations(tmp_out="./tmp",dyn_yamlfile=yamlfile_e,cdynafile=cdyna_e,teth5file=teth5_e,ndbfile=save_path+'/SAVE/'+ndb)
# read hole
hole_occups = dynamic_occupations(tmp_out="./tmp",dyn_yamlfile=yamlfile_h,cdynafile=cdyna_h,teth5file=teth5_h,ndbfile=save_path+'/SAVE/'+ndb)
elec_occups.pert_grid_reduced()
hole_occups.pert_grid_reduced()


ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')
y_k_grid=ylat.k_grid
p_k_grid=elec_occups.kpts_grid

### Generate fake grids to check the code ##############
if Debug:
    y_grid=np.array([4,4,1])
    y_k_grid=y_grid

p_grid=np.array([24,24,24])
p_k_grid=p_grid

#######################################################


print("Perturbo k-grid ",p_k_grid)
print("Yambo k-grid ",y_k_grid)

if any(p_k_grid%y_k_grid != 0):
    print("Incompatible k-grid!! ")
    sys.exit(0)
else:
    print("Compatible k-grids found")
    
# This code works only for ratio%2
grid_ratio=np.int32(p_k_grid/y_k_grid)
print("Grid ratio ",grid_ratio)
if(any(grid_ratio%2 ==0)):
   print("Please use an odd grid ratio! ")
   sys.exit(0)

elec_occups.read_pert_num_kpts()
if Debug:
    n_kpt_pert=np.prod(p_k_grid)
else:
    n_kpt_pert=elec_occups.num_pert_kpts
print("Number of k-points in perturbo : ",n_kpt_pert)

if Debug:
    elec_pert_kpts    =generate_grid(p_k_grid)
    hole_pert_kpts    =generate_grid(p_k_grid)
    yambo_kpts_bz=generate_grid(y_k_grid)
    yambo_kpts_ibz=generate_grid(y_k_grid)
else:
    elec_pert_kpts    =elec_occups.read_perturbo_kpts()
    hole_pert_kpts    =hole_occups.read_perturbo_kpts()
#   elec_pert_kpts    =generate_grid(p_k_grid)
#   hole_pert_kpts    =generate_grid(p_k_grid)
    yambo_kpts_bz =ylat.red_kpoints
    yambo_kpts_ibz=ylat.get_ibz_kpoints(units='red')


# Bring all k-points between [0,1)
elec_pert_kpts=make_kpositive(elec_pert_kpts.tolist())
hole_pert_kpts=make_kpositive(hole_pert_kpts.tolist())
yambo_kpts_ibz=make_kpositive(yambo_kpts_ibz.tolist())

# Find the smallest q-vector in crystal coordinates
small_q=np.zeros(3,dtype=float)
small_q=1.0/p_k_grid

# Build integer k-points for perturbo (fine grid)
elec_pert_ikpt=[np.int32(np.rint(kpt/small_q)) for kpt in elec_pert_kpts]
hole_pert_ikpt=[np.int32(np.rint(kpt/small_q)) for kpt in hole_pert_kpts]
with open('perturbo_ik_bz.pts', 'w') as f:
    f.write("#Perturbo elec ik-points in the BZ\n")
    for ikpt in elec_pert_ikpt:
        f.write(str(ikpt)+'\n')
    f.write("#Perturbo hole ik-points in the BZ\n")
    for ikpt in hole_pert_ikpt:
        f.write(str(ikpt)+'\n')

# Bulkd integer BZ k-points for Yambo grid using the q-vector of perturbo
yambo_ikpt=[np.int32(np.rint(kpt/small_q)) for kpt in yambo_kpts_bz]
with open('yambo_ik_bz.pts', 'w') as f:
    f.write("#Yambo ik-points in the BZ\n")
    for ikpt in yambo_ikpt:
        f.write(str(ikpt)+'\n')

print("Number of IBZ k-points in Yabmo : ",len(yambo_kpts_ibz))

# Bulkd integer IBZ k-points for Yambo (course grid)
yambo_ikpt_ibz=[np.int32(np.rint(kpt/small_q)) for kpt in yambo_kpts_ibz]
with open('yambo_ik_ibz.pts', 'w') as f:
    f.write("#Yambo ik-points in the IBZ\n")
    for ikpt in yambo_ikpt_ibz:
        f.write(str(ikpt)+'\n')

#search perturbo neighboars for each yambo point in the IBZ
print("Find electron neighboars :")
elec_yneighboars=find_neighboars(yambo_ikpt_ibz,elec_pert_ikpt,p_k_grid,grid_ratio,Debug)
print("Find hole neighboars :")
hole_yneighboars=find_neighboars(yambo_ikpt_ibz,hole_pert_ikpt,p_k_grid,grid_ratio,Debug)

# Get bands indexes
elec_occups.get_vcb_indices()
elec_occups.parse_bands_from_yaml()
hole_occups.get_vcb_indices()
hole_occups.parse_bands_from_yaml()

print("Electrons bands : ",elec_occups.pert_bands)
print("Hole bands : ",hole_occups.pert_bands)

# Read occupations
hole_occups.get_files()
elec_occups.get_files()

# Copy bare occupation for all k-points
# read Yambo DB
RT_db=YamboRT_Carriers_DB(calc=save_path+'/SAVE/',carriers_db='ndb.RT_carriers')
RT_db.get_info()

# Make a folder to store carrries files
carriers_path="CARRIERS" 
if not os.path.exists(carriers_path):
    os.makedirs(carriers_path)

#Array of the occupation variation
delta_f=np.zeros([elec_occups.num_kpts,elec_occups.full_num_bands],float)

#Copy originak ndb.RT_carries in the folder and then modify it
for occ_key in elec_occups.occupation.keys():
    shutil.copyfile(save_path+'SAVE/ndb.RT_carriers',carriers_path+"/ndb.RT_carriers_"+str(occ_key))
    # Open database without closing it "keep_open=True"
    RT_db_new=YamboRT_Carriers_DB(calc=carriers_path,carriers_db='ndb.RT_carriers_'+str(occ_key),keep_open=True)

    copy_occupation(elec_yneighboars,elec_occups,occ_key,delta_f,'e')
    copy_occupation(hole_yneighboars,hole_occups,occ_key,delta_f,'h')
    
    # Update data in the database
    RT_db_new.delta_f=np.reshape(delta_f,elec_occups.num_kpts*elec_occups.full_num_bands)
    RT_db_new.updateDB()

    # Close the database
    RT_db_new.closeDB()
    
