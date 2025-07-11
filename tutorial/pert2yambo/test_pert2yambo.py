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



dynoccups = dynamic_occupations(tmp_out="./tmp",dyn_yamlfile=yamlfile_e,cdynafile=cdyna_e,teth5file=teth5_e,ndbfile=save_path+'/SAVE/'+ndb)
dynoccups.pert_grid_reduced()

ylat = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1')
y_k_grid=ylat.k_grid
p_k_grid=dynoccups.kpts_grid

print("Perturbo k-grid ",p_k_grid)
print("Yambo k-grid ",y_k_grid)

if any(p_k_grid%y_k_grid != 0):
    print("Incompatible k-grid!! ")
    sys.exit(0)
else:
    print("Compatible k-grids found")

dynoccups.read_pert_num_kpts()
n_kpt_pert=dynoccups.num_pert_kpts
print("Number of k-points in perturbo : ",n_kpt_pert)



pert_kpts=dynoccups.read_perturbo_kpts()

pert_kpts=make_kpositive(pert_kpts.tolist())
small_q=np.zeros(3,dtype=float)
small_q=1.0/p_k_grid

pert_ikpt=[np.int32(np.rint(kpt/small_q)) for kpt in pert_kpts]
# for ikpt in pert_ikpt:
#    print(ikpt)


yambo_ikpt=[np.int32(np.rint(kpt/small_q)) for kpt in ylat.red_kpoints]
for ikpt in yambo_ikpt:
    print(ikpt)
# dynoccups.get_files()


dynoccups.get_vcb_indices()
# dynoccups.parse_bands_from_yaml()
# dynoccups.get_files()

