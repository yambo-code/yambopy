from yambopy import *
from pylab import *

lattice = YamboLatticeDB(save='rt/SAVE')

fig = plt.figure()
a = YamboElectronPhononDB(lattice,filename='ndb.elph_gkkp_expanded',folder_gkkp='rt/GKKP',save='rt/SAVE')
a.plot_map(fig,ib1=1,ib2=1,ik1=7,all_phonons=False,size=20,lim=0.25)
plt.show()

fig = plt.figure()
a.plot_map(fig,ib1=1,ib2=1,ik1=7,all_phonons=True,size=100,lim=0.40)
plt.show()
