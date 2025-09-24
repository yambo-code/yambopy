from qepy import *

uc="diamond.scf.in"
eigv_filename ="matdyn.modes"
qe_input = PwIn.from_file(uc)
qe_dyn   = Matdyn.from_modes_file(filename=eigv_filename)

R=[2,2,2]

Map_Phonons(qe_input, qe_dyn, R)
