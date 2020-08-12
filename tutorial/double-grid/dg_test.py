from qepy import *
from yambopy import *
from materials import *

# Path of executables
pw = 'pw.x'
p2y = '/Users/fulvio.paleari/software/yambo-andreaM/bin/p2y'
yambo = '/Users/fulvio.paleari/software/yambo-andreaM/bin/yambo'
ypp = '/Users/fulvio.paleari/software/yambo-andreaM/bin/ypp'

prefix = 'bn'

# List of coarse grids (CG)
cg_grids = [[3,3,1],[6,6,1],[9,9,1],[12,12,1]]
# List of random fine grids (FG)
fg_grids = [['9_fg','18_fg','36_fg'],['36_fg','72_fg','144_fg'],['81_fg'],['144_fg']]

# Laser energy (eV)
E_laser = 5.

# Paths of input data: QE scf save, pseudos
input_data = '/Users/fulvio.paleari/software/whypy/yambo-whypy-devel/tutorial/double-grid/input_data'
scf_save_path = input_data
pseudo_path = '%s/pseudos'%input_data
work_dir = '/Users/fulvio.paleari/software/whypy/yambo-whypy-devel/tutorial/double-grid'

#Submission script to scheduler
qe_run_script = None
yambo_run_script = None

# BN inputs
bn_inp = hBN_1l_test(prefix=prefix,pseudo_path=pseudo_path)
qe_input    = bn_inp.nscfin
yambo_input = bn_inp.ipin

#QE output prefix
nscf_out = "nscf" #"slurm"
#Yambo output folder(s):
y_out_dir = "results"

#Call to main class
YamboDG_Optimize(cg_grids,fg_grids,prefix,qe_input,yambo_input,scf_save_path,pseudo_path,RUN_path=work_dir,nscf_out=nscf_out,y_out_dir=y_out_dir,E_laser=E_laser,pw=pw,yambo=yambo,ypp=ypp,p2y=p2y,STEPS='all')
