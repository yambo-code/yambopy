# This script reads the carrier database
# and display it along a path in histogram form

from yambopy import *
import matplotlib.gridspec as gridspec

#save = 'rt-24x24/SAVE' # Save with the appropriate band structure
folder = 'rt-24x24'
calc   = 'QSSIN-D-100.0fs-2.07eV-300K-DG' # Where RT carrier output is
source = '%s/%s/pulse/'%(folder,calc)
path = [[0.0,0.0,0.0],[0.5,0.0,0.0],[0.33333,0.33333,0.0],[0.0,0.0,0.0]]
nbv = 2 ; nbc = 2 # nb of valence and conduction bands
degen_thresh = 0.050 # from what diff in E do we consider degen states


# Instances containing LDA bandstructure and occupations
yrt = YamboRTDB(folder=folder,calc=calc)
yrt.get_path(path)

# aliases
kindex = yrt.bands_indexes # kpoint indexes (in order) to draw path
eigenvalues = yrt.eigenvalues[kindex,:] # eigenvalues (LDA) of the band included in the RT simulation

nbands = yrt.nbands # number of bands in the RT simulation
times = [i * 1e15 for i in yrt.times] # times are now in fs
occupations = yrt.occupations[:,kindex,:] # format time,kindex,band index (from 0 to nbands)

if nbv+nbc != nbands:
    raise NameError('Incompatible number of bands, set nbv and nbc in script.')



# occupations in CBs/VBs are summed for easier reading if there are more than one
if nbv > 1:
    # one entry per band +1 for the total occ
    occ_v = np.zeros((len(times),len(kindex),nbv+1))
    occ_c = np.zeros((len(times),len(kindex),nbc+1))
    for n in range(nbv):
        occ_v[:,:,n] = -occupations[:,:,n] # minus sign to get positive occupations
        np.add(occ_v[:,:,n],occ_v[:,:,nbv],occ_v[:,:,nbv]) # each time we add the occ of the current band to the total
    for n in range(nbc):
        occ_c[:,:,n] = occupations[:,:,n+nbv] # +nbv to read CBs
        np.add(occ_c[:,:,n],occ_c[:,:,nbv],occ_c[:,:,nbc]) # each time we add the occ of the current band to the total


#    for k in kindex:
#        if abs(occupations[0,k,0]-occupations[0,k,1])<degen_thres:
#            occ_v = occupations[:,:,0]+occupations[:,:,1]

# The external field is read from the o- file
ext = np.loadtxt(source+'o-pulse.external_field')
field = ext[:,2]/max(abs(ext[:,2])) # polarization : x=1,y=2,z=3

# Gridspec allows to place subplots on a grid
# spacing for exemple can be customised
gs = gridspec.GridSpec(8, 1)
gs.update(hspace=1.1) # bigger horizontal spacing

# To plot occ properly, we need an order array of integers
xocc = [i for i in range(len(kindex))]


#for t in range(len(times)):
for t in (10,):

    ax1 = plt.subplot(gs[:-1, :]) # bandstructure w/ occupation plot
    plt.xticks(()) # remove x ticks
# Band structure
    plt.plot(eigenvalues,'k-')
# occupation in the form of histogram
    plt.bar(xocc,occ_v[t,:,nbv]*1000,width=0.02,bottom=eigenvalues[:,nbv-1],color='blue',edgecolor='none')
    plt.bar(xocc,occ_c[t,:,nbc]*1000,width=0.02,bottom=eigenvalues[:,nbands-1],color='red',edgecolor='none')

    plt.ylabel('E (eV)')
    plt.text(0.50,0.5, '%d fs'%times[t])


    ax2 = plt.subplot(gs[-1, :])
    plt.xticks(())
    plt.yticks(())
    plt.ylabel('Field')
    plt.xlim((0,times[-1]))
    plt.ylim((-1.3,1.3))
    plt.plot(field[:int(times[t])])
    plt.show()

