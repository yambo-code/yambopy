#
# This example show how to calculate Density of States using YamboPy
#
from qepy import *
from yambopy import *
import matplotlib.pyplot as plt

save_folder='./SAVE' 
yel = YamboElectronsDB.from_db_file(folder=save_folder)
dos=yel.getDOS()

fig = plt.figure(figsize=(6,4))
ax.set_xlabel('$Energy[eV]$')
ax.set_ylabel('DOS$')
plt.
