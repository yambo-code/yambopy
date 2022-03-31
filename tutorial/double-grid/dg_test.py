from qepy import *
from yambopy import *
from materials import *

# Path of executables
pw_path = '/Users/fulvio.paleari/software/q-e/bin'
y_path = '/Users/fulvio.paleari/software/yambo-andreaM/bin'

prefix = 'bn'

# List of coarse grids (CG)
#cg_grids = [[3,3,1],[6,6,1],[9,9,1],[12,12,1]]
cg_grids = [[3,3,1],[6,6,1]]
# List of random fine grids (FG)
#fg_grids = [[9,18,36],[36,72,144],[81],[144]]
fg_grids = [[9,18,36],[36]]

# Laser energy (eV)
E_laser = 5.

# Paths of input data: QE scf save, pseudos
input_data = '/Users/fulvio.paleari/software/whypy/yambo-whypy-devel/tutorial/double-grid/input_data'
scf_save_path = input_data
pseudo_path = '%s/pseudos'%input_data
work_dir = '/Users/fulvio.paleari/software/whypy/yambo-whypy-devel/tutorial/double-grid'

# BN inputs
bn_inp = hBN_1l_test(prefix=prefix,pseudo_path=pseudo_path)
qe_input    = bn_inp.nscfin
yambo_input = bn_inp.ipin

#QE output prefix
nscf_out = "nscf" #"slurm"
#Yambo output folder(s):
y_out_dir = "results"

#Call to main class [convergence on]
YamboDG_Optimize(cg_grids,fg_grids,prefix,qe_input,yambo_input,E_laser=E_laser,STEPS='all',converge_DG=True,\
                scf_save_path=scf_save_path,pseudo_path=pseudo_path,RUN_path=work_dir,\
                nscf_out=nscf_out,y_out_dir=y_out_dir,\
                pw_exec_path=pw_path,yambo_exec_path=y_path,yambo_exec='yambo',\
                save_type='simple',yambo_calc_type='ip')
