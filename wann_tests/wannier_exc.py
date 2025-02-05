import tbmodels
from yambopy import *
import matplotlib.pyplot as plt
from yambopy.lattice import car_red, red_car
import matplotlib.pylab as pylab
from pylab import rcParams
from matplotlib.ticker import MultipleLocator
import matplotlib.lines as mlines
# Set figure size and DPI globally
rcParams['figure.figsize'] = 6, 4  # 6 inches wide, 4 inches tall
rcParams['figure.dpi'] = 300  # 300 dots per inch

# Font sizes
rcParams['font.size'] = 20
rcParams['axes.titlesize'] = 20  # fontsize of the axes tit bels
rcParams['xtick.labelsize'] = 20  # fontsize of the tick labels
rcParams['ytick.labelsize'] = 20  # fontsize of the tick labels
rcParams['legend.fontsize'] = 16  # legend fontsize

# Line widths and marker sizes
rcParams['lines.linewidth'] = 2  # width of lines
rcParams['lines.markersize'] = 6  # size of markers

# Label, axis, and tick thickness
rcParams['axes.linewidth'] = 1.5  # axis line width
rcParams['xtick.major.width'] = 1.2  # width of the major tick lines
rcParams['ytick.major.width'] = 1.2  # width of the major tick lines

WORK_PATH='/gpfs/work2/0/prjs0229/share/sb-rr/LiF/July2024/'
QE_PATH='/gpfs/work2/0/prjs0229/share/sb-rr/LiF/July2024/myqe/'
YAMBO_PATH = f'/gpfs/work2/0/prjs0229/share/sb-rr/LiF/July2024/Optics/yambo-DS/'
ry2ev=13.605703976

# Istance of useful classes
savedb_k = YamboSaveDB.from_db_file(f'{YAMBO_PATH}/12x12x12/SAVE')
lat_k = YamboLatticeDB.from_db_file(f'{YAMBO_PATH}/12x12x12/SAVE/ns.db1')
savedb_q = YamboSaveDB.from_db_file(f'{YAMBO_PATH}/6x6x6/SAVE')
lat_q = YamboLatticeDB.from_db_file(f'{YAMBO_PATH}/6x6x6/SAVE/ns.db1')

# Create instance of real space Hamiltonian in MLWF basis
hrk=HR(f'{QE_PATH}/nscf-wannier-12x12x12/LiF')
hrq = HR(f'{QE_PATH}/nscf-wannier-6x6x6/LiF')
#hrqexc=HR(f'{YAMBO_TUT_PATH}/unshifted-grid/nscf-wannier-kgrid/exc/LiF_exc')

nnkp_kgrid = NNKP_Grids(f'{QE_PATH}/nscf-wannier-12x12x12//LiF', lat_k, yambo_grid=True)
nnkp_qgrid = NNKP_Grids(f'{QE_PATH}/nscf-wannier-6x6x6//LiF', lat_q, yambo_grid=True)
print('start grids')
# We need all these auxiliary grids for wannierization of the BSE Hamiltonian
nnkp_kgrid.get_kmq_grid(nnkp_qgrid)
nnkp_kgrid.get_qpb_grid(nnkp_qgrid)
nnkp_qgrid.get_qpb_grid(nnkp_qgrid)
nnkp_kgrid.get_kpbover2_grid(nnkp_qgrid)
nnkp_kgrid.get_kmqmbover2_grid(nnkp_qgrid)
print('end grids')
model = TBMODEL.from_wannier_files(
    hr_file=f'{QE_PATH}/nscf-wannier-12x12x12/LiF_hr.dat',
    wsvec_file=f'{QE_PATH}/nscf-wannier-12x12x12/LiF_wsvec.dat',
    #xyz_file=f'{WORK_PATH}/w90files/hBN_centres.xyz',
    win_file=f'{QE_PATH}/nscf-wannier-12x12x12/LiF.win'
)
model.set_mpgrid(nnkp_kgrid)

fermie = 1.0
model.solve_ham_from_hr(lat_k, hrk, fermie=fermie  )

print ('before h2p')
h2p = H2P(model, f'{YAMBO_PATH}/12x12x12/SAVE',qmpgrid=nnkp_qgrid,bse_nv=3,bse_nc=1,
          cpot=None,kernel_path =f'{YAMBO_PATH}/12x12x12/06_BSE_triplet/',
          excitons_path = f'{YAMBO_PATH}/12x12x12/06_BSE_triplet/',TD=False, method='kernel',run_parallel=True)

print('built h2p')
h2p.solve_H2P()

print('writing files')

h2p.write_exc_overlap(seedname='LiF_exc',trange=np.arange(0,3), tprange=np.arange(0,3))
h2p.write_exc_eig(seedname='LiF_exc',trange=np.arange(0,3))
h2p.write_exc_nnkp(seedname='LiF_exc',trange=np.arange(0,3))
h2p.write_exc_amn(infile=f'{YAMBO_PATH}/May2024/nscf-wannier-12x12x12/LiF',seedname='LiF_exc',trange=np.arange(0,3), tprange=np.arange(0,3))

