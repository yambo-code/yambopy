#
# This example show how to calculate Density of States using YamboPy
#
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt

save_folder='./SAVE' 
yel = YamboElectronsDB.from_db_file(folder=save_folder)
energies,dos=yel.getDOS(emin=-14.0,emax=20.0)

fig = plt.figure(figsize=(6,4))
plt.xlabel('$Energy[eV]$')
plt.ylabel('DOS$')
plt.plot(energies,dos)
plt.show()
