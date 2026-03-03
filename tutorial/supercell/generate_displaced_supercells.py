from supercell_utilities import *
from yambopy import *


qe_dyn = Matdyn.from_modes_file(filename='qe/matdyn.modes_full_3x3')
qe = PwIn.from_file('qe/AlN.scf.in')
R = [3, 3, 1]

sc = MySupercell(qe)

sc.apply_full_thermal_displacement(R, qe_dyn, T_Kelvin=0, rm_output_sym=True, 
                                   qe_control_dict={'pseudo_dir':"'../pseudo/'"}, 
                                   qe_system_dict={'nbnd':80},
                                   output_all=True
                                   )

