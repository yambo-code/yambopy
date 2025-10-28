from qepy import *

uc="diamond.scf.in"
eigv_filename ="matdyn.modes"
qe = PwIn.from_file(uc)
qe_dyn = Matdyn.from_modes_file(filename=eigv_filename)

generate_ZG_conf(qe, qe_dyn, modes=[3,4,5])
