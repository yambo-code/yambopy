from qepy import *

uc="diamond.scf.in"
eigv_filename ="matdyn.modes"
qe = PwIn.from_file(uc)
qe_dyn = Matdyn.from_modes_file(filename=eigv_filename)

