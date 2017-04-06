# This script reads the carrier database
# and display it along a path in histogram form

from yambopy import *
import matplotlib.gridspec as gridspec

save = 'rt-24x24/SAVE' # Save with the appropriate band structure
source = 'rt-24x24/QSSIN-100.0fs-2.07eV-300K/pulse/' # Where RT carrier output is
path = [[0.0,0.0,0.0],[0.5,0.0,0.0],[0.33333,0.33333,0.0],[0.0,0.0,0.0]]


# Instances containing LDA bandstructure and occupations
ys = YamboSaveDB(save)
yrt = YamboRTDB(ys,path=source)
ys.get_path(path)

# aliases
kindex = ys.bands_indexes # kpoint indexes (in order) to draw path
eigenvalues = ys.eigenvalues[kindex,yrt.nband_min-1:yrt.nband_max] # eigenvalues (LDA) of the band included in the RT simulation
nbands = yrt.nbands # number of bands in the RT simulation
times = [i * 1e15 for i in yrt.times] # times are now in fs
occupations = yrt.occupations[:,kindex,:] # format time,kindex,band index (from 0 to nbands)

# The external field is read from the o- file
ext = np.loadtxt(source+'o-pulse.external_field')
field = ext[:,2]/max(abs(ext[:,2])) # polarization : x=1,y=2,z=3

# Gridspec allows to place subplots on a grid
# spacing for exemple can be customised
gs = gridspec.GridSpec(8, 1)
gs.update(hspace=1.1) # bigger horizontal spacing

# To plot occ properly, we need an order array of integers
xocc = [i for i in range(len(kindex))]


# plot at time 10 for exemple
t = 10

ax1 = plt.subplot(gs[:-1, :]) # bandstructure w/ occupation plot
plt.xticks(()) # remove x ticks
# Band structure
plt.plot(eigenvalues,'k-')
# occupation in the form of histogram
nb = 0
plt.bar(xocc,(-1)**nb*occupations[t,:,nb]*1000,width=0.02,bottom=eigenvalues[:,nb],color=['blue'],edgecolor='none')
nb = 1
plt.bar(xocc,(-1)**nb*occupations[t,:,nb]*1000,width=0.02,bottom=eigenvalues[:,nb],color='blue',edgecolor='none')
nb=2
plt.bar([i for i in range(len(kindex))],(-1)**(nb+1)*occupations[t,:,nb]*1000,width=0.02,bottom=eigenvalues[:,nb],color='red')
nb=3
plt.bar([i for i in range(len(kindex))],(-1)**(nb+1)*occupations[t,:,nb]*1000,width=0.02,bottom=eigenvalues[:,nb],color="red")

plt.ylabel('E (eV)')
plt.text(0.50,0.5, '%d fs'%times[t])


ax2 = plt.subplot(gs[-1, :])
plt.xticks(())
plt.yticks(())
plt.ylabel('Field')
plt.xlim((0,times[-1]))
plt.ylim((-1.2,1.2))
plt.plot(field[:int(times[t])])
plt.show()
exit()


# plot at time 10 for exemple for band "0"
t = 10
nb = 0

plt.show()
exit()

# we sti
plt.plot(ys.eigenvalues[ys.bands_indexes,yrt.nband_min-1:yrt.nband_max],'k-')
plt.bar([i for i in range(len(ys.bands_indexes))],[i for i in range(len(ys.bands_indexes))],width=0.02,bottom=ys.eigenvalues[ys.bands_indexes,yrt.nband_min-1])
plt.show()
