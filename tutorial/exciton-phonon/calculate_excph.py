import numpy as np
from yambopy import YamboLatticeDB,YamboWFDB,LetzElphElectronPhononDB
from yambopy.exciton_phonon.excph_matrix_elements import exciton_phonon_matelem

#path = '3D_hBN'
#bands_range=[6,10]
path = '1L_MoS2'
bands_range=[24,28] # 2 valence bands, 2 conduction bands
nexc = 12 # number of excitonic states

bsepath    = f'{path}/bse-allq_full' # Path to BSE calculation (Lin=Lout)
savepath   = f'{path}/SAVE'     # Yambo SAVE
ndb_elph   = f'{path}/ndb.elph' # LetzElPhC electron-phonon database (any convention)

# Read lattice
lattice = YamboLatticeDB.from_db_file(filename=f'{savepath}/ns.db1')
# Read electron-phonon
elph    = LetzElphElectronPhononDB(ndb_elph,read_all=False)
# Read wave functions
wfcs    = YamboWFDB(filename='ns.wf',save=savepath,latdb=lattice,bands_range=bands_range)
# Calculate exciton-phonon matrix elements
exph = exciton_phonon_matelem(lattice,elph,wfcs,BSE_dir=bsepath,neigs=nexc,dmat_mode='save',exph_file='MoS2_Ex-ph.npy')
