from yambopy import *
import tbmodels
import matplotlib.pyplot as plt
from yambopy.lattice import car_red, red_car
import matplotlib.pylab as pylab
from pylab import rcParams
from matplotlib.ticker import MultipleLocator
import matplotlib.lines as mlines


shared_dir = '/gpfs/work2/0/prjs0229/share/sb-rr/WS2/'
savedb_k = YamboSaveDB.from_db_file(f'{shared_dir}nscf-6x6x1/out/ws2.save/SAVE')
lat_k = YamboLatticeDB.from_db_file(f'{shared_dir}nscf-6x6x1/out/ws2.save/SAVE/ns.db1')
hrk=HR(f'{shared_dir}/nscf-wannier-6x6x1/ws2')
hrq=HR(f'{shared_dir}/nscf-wannier-3x3x1/ws2')

nnkp_kgrid = NNKP_Grids(f'{shared_dir}/nscf-wannier-6x6x1/ws2', lat_k, yambo_grid=True)
nnkp_qgrid = NNKP_Grids(f'{shared_dir}/nscf-wannier-3x3x1/ws2', lat_k, yambo_grid=True)

# # We need all these auxiliary grids for wannierization of the BSE Hamiltonian
nnkp_kgrid.get_kmq_grid(nnkp_qgrid)
nnkp_kgrid.get_qpb_grid(nnkp_qgrid)
nnkp_qgrid.get_qpb_grid(nnkp_qgrid)
nnkp_kgrid.get_kpbover2_grid(nnkp_qgrid)
nnkp_kgrid.get_kmqmbover2_grid(nnkp_qgrid)

model = TBMODEL.from_wannier_files(
    hr_file=f'{shared_dir}nscf-wannier-6x6x1/ws2_hr.dat',
    wsvec_file=f'{shared_dir}nscf-wannier-6x6x1/ws2_wsvec.dat',
    #xyz_file=f'{WORK_PATH}/w90files/hBN_centres.xyz',
    win_file=f'{shared_dir}nscf-wannier-6x6x1/ws2.win'
)
#set the grid
model.set_mpgrid(nnkp_kgrid)


npoints = 30
path_kpoints = Path([[[  0.0,  0.0,  0.0],'$\Gamma$'],
                     [[0.5, 0.0,  0.0],'M'],  
              [[  1/3,  1/3,  0.0],'K'],
              [[  0.0,  0.0,  0.0],'$\Gamma$']],[npoints,npoints,npoints] )
klist = path_kpoints.get_klist()

kpoints_red = path_kpoints.get_klist()[:,0:3]
kpoints_car = red_car(kpoints_red, lat_k.rlat)
kdistance = path_kpoints.distances
kpoints = path_kpoints.kpoints
H_atk = model.solve_ham_from_hr(lat_k, hrk , fermie = 1.0)
E_k_tb= np.array(model.eigenval(kpoints_red))
E_k= np.array(model.get_eigenval(kpoints_red,from_hr=True))


fermie = -4.0
model.solve_ham_from_hr(lat_k, hrk, fermie=fermie)
hlm = model.get_hlm(lat_k.lat, hrk)





savedb_q = YamboSaveDB.from_db_file(f'{shared_dir}/nscf-wannier-3x3x1/out/ws2.save/SAVE')
lat_q = YamboLatticeDB.from_db_file(f'{shared_dir}/nscf-wannier-3x3x1/out/ws2.save/SAVE/ns.db1')

# fermie = -4
# model.solve_ham_from_hr(lat_k, hrk, fermie=fermie)

cpot = CoulombPotentials(ngrid=nnkp_kgrid, lattice=lat_k, v0=1.0, lc=15.0, w=1.0, r0=1.0, tolr=0.001, ediel=[1.0,1.0,1.0])

savedb_path = f'{shared_dir}/nscf-wannier-6x6x1/out/ws2.save/SAVE'
kernel_path = f'{shared_dir}bse/out/ws2.save/ws2BSE/'
h2p = H2P(model, savedb_path, qmpgrid=nnkp_kgrid, kernel_path=kernel_path, method='model',ctype='v2drk', cpot=cpot)

h2p.solve_H2P()

w, eps = h2p.get_eps(hlm=model.hlm, emin=0.0, emax=10, estep=0.1, eta=1., with_bse=False)


print('Writing files')

np.save('eps.npy', np.array([w, eps[:,0,1]]))