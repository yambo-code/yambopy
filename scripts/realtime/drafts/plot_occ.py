# This script reads the carrier database
# and display it along a path in histogram form

from yambopy import *
import matplotlib.gridspec as gridspec


############
# SETTINGS #
############

#save = 'rt-24x24/SAVE' # Save with the appropriate band structure
folder = 'rt-24x24'
#calc   = 'QSSIN-D-100.0fs-2.07eV-300K-DG' # Where RT carrier output is
calc   = 'QSSIN-100.0fs-2.07eV-300K' # Where RT carrier output is
source = '%s/%s/pulse/'%(folder,calc)
path = [[0.0,0.0,0.0],[0.5,0.0,0.0],[0.33333,0.33333,0.0],[0.0,0.0,0.0]]
nbv = 2 ; nbc = 2 # nb of valence and conduction bands
occ_scaling = 1 # max occupation will be 1eV high


########
# INIT #
########

# Instances containing LDA bandstructure and occupations
yrt = YamboRTDB(folder=folder,calc=calc)
yrt.get_path(path)

# aliases
kindex = yrt.bands_indexes # kpoint indexes (in order) to draw path
eigenvalues = yrt.eigenvalues[kindex,:] # eigenvalues (LDA) of the band included in the RT simulation

nbands = yrt.nbands # number of bands in the RT simulation
times = [i * 1e15 for i in yrt.times] # times are now in fs
occupations = yrt.occupations[:,kindex,:]/np.amax(yrt.occupations[:,kindex,:])*occ_scaling # format time,kindex,band index (from 0 to nbands)

if nbv+nbc != nbands:
    raise NameError('Incompatible number of bands, set nbv and nbc in script.')

##################
# ENERGY DISTRIB #
##################

### Idea
## Now that histogram is mastered, I can split eigenvalues and occupations in CB and VB before doing the histogram, so that I can automatically use bins=len(array)

# yrt.occupations[t,k,n]
# yrt.eigenvalues[k,n]
i=0;j=0
list_e=[] ; list_h=[]
for k in range(yrt.nkpoints):
    for n in range(yrt.nbands):
        e = yrt.eigenvalues[k,n]
        if e<=0.0:
            list_h.append((k,n))
        else:
            list_e.append((k,n))

occ_e = np.zeros((len(times),len(list_e),2))
for t in range(len(times)):
    for i,(k,n) in enumerate(list_e):
        occ_e[t,i,0]=yrt.eigenvalues[k,n]
        occ_e[t,i,1]=yrt.occupations[t,k,n]

occ_h = np.zeros((len(times),len(list_h),2))
for t in range(len(times)):
    for i,(k,n) in enumerate(list_h):
        occ_h[t,i,0]=yrt.eigenvalues[k,n]
        occ_h[t,i,1]=yrt.occupations[t,k,n]

####################################################
###### Comparison between ypp and python data ######
#ypp = np.loadtxt('rt-24x24/QSSIN-100.0fs-2.07eV-300K/o-pulse.YPP-RT_occupations_DATA_1_of_7').T
#plt.scatter(ypp[0,:],ypp[10,:],color='b',marker='o')        # The columns value have the same data
#plt.scatter(occ_e[9,:,0],occ_e[9,:,1],color='r',marker='+') # So be careful w/ index if willing to do comparison
#plt.show()
#
#plt.hist(occ_e[9,:,0], weights=occ_e[9,:,1],bins=len(occ_e[0]),color='b',edgecolor='none') # Not sure if ypp sums the almost degen C bands at
#plt.hist(ypp[0,:]+0.1,weights=ypp[10,:],bins=len(ypp[0,:]),color='r',edgecolor='none')     # K or if it is artifact of plotting
#plt.show()
####################################################

#### Latest question : should almost degenerate states be merged for the fitting ? What would be the energy threshold ?
### Answer : no, use Counter to sum exactly (thr 10e-8) degenerate states and then plot normally

# occupation as f(E). occe[t] is (e,occ(t))
#occe = [(yrt.eigenvalues.flatten(),occt.flatten()) for occt in yrt.occupations]
## tuple (e,occ) @ t is (occe[t][0][i],occe[t][1][i])
#occe =  np.array(occe)

## Removing degeneracies
## case for t 10
#dic = {}
#for i in range(len(occe[10][0])):
#    if occe[10][0][i] not in dic:
#        dic[ occe[10][0][i] ] = 0
#    dic[ occe[10][0][i] ] += occe[10][1][i]
#
#result = list(dic.items())
#print len(result)
#print result[0]
#
#plt.scatter(result)
#plt.show()
#exit()
#
##occetest = occe[10]
#print occetest[0][0],occetest[1][0]
#print occetest[0][1],occetest[1][1]
#print occetest[0][2],occetest[1][2]
#print len(occetest)
#print len(occetest[0])
##occetest = sorted(occetest)
##print occetest
#plt.scatter(occetest[0],occetest[1])
#plt.show()
#
#sum = {}
#print len(occetest)
#i = 0
#for eocc in occetest:
#    i+=1
#    if not eocc[0] in sum:
#        sum[eocc[0]]=0
#    sum[eocc[0]]+=eocc[1]
#print i
#print sum
#result = list(sum.items())
#print result
#plt.scatter(result)
#plt.show()
#exit()
##np.savetxt('occe.dat',occe[10])
#print occe
#exit()
###########
## I can use scatter to plot, but still I have not solved the degeneracy problem
## Either I put energies in order to have a correct histogram to which I can pass the energies in order to make the right bins
## => Or I sum degenerate energies' occupations and I can scatter and fit easily
#for t in range(len(times)):
#    plt.scatter(occe[t][0],occe[t][1])
#    plt.show()
t=10
print occe[t].shape
print occe[t].T.shape
#testarray = np.where( occe[t][0,:] <= 0.0)
#print testarray
#plt.hist(testarray[0],weights=testarray[1],bins=len(testarray[0]))
#plt.show()
hist = np.histogram(occe[t][0,:],weights=occe[t][1,:],bins=1200)
print hist[0].shape
#plt.hist(occe[t].T[:,0], weights=occe[t].T[:,1],bins=1200)  # plt.hist passes it's arguments to np.histogram
plt.hist(occe[t][0,:], weights=occe[t][1,:],bins=1200)  # plt.hist passes it's arguments to np.histogram
plt.show()
exit()

#occe = np.zeros((yrt.nkpoints*nbands)) # occ_energy
#occo = np.zeros((len(times),yrt.nkpoints*nbands)) # occ_occ
#
#occe = yrt.eigenvalues.flatten()
#for t in range(len(times)):
#    occo[t]=yrt.occupations[t].flatten()
#
#array = np.stack((occo,occe),axis=-1)
#print array
#exit()
#plt.hist(occo[10], bins=occe)  # plt.hist passes it's arguments to np.histogram
#plt.show()
#exit()

## I want an array with, for each time, the occupation at a given energy level
## The energies are constants
## First, find (k,n) combinations that give degenerate energy levels
## Degenerate eps => two id values in yrt.eigenvalues
#degen_eigen = []
#it = np.nditer(yrt.eigenvalues, flags=['multi_index']) # numpy way to go through an array
#while not it.finished:
#    # it[0] is current value
#    # it.multi_index is tuple (k,n)
#    if any(yrt.eigenvalues[it.multi_index] in tup for tup in degen_eigen):
#        # already exists, just expand current tuple
#        tup.append(it.multi_index)
#    else:
#        degen_eigen.append((yrt.eigenvalues[it.multi_index],it.multi_index))
#    it.iternext() # next tuple (k,n)
## Then, build the array
print degen_eigen
exit()

#################
# BAR PLOT DATA #
#################

# occupations in CBs/VBs are summed for easier reading if there are more than one
#TODO add a color scheme to see what band contributes the most to the occupation at a k point
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


##############
# EXT. FIELD #
##############

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
    plt.bar(xocc,occ_v[t,:,nbv],width=0.5,bottom=eigenvalues[:,nbv-1],color='blue',edgecolor='none')
    plt.bar(xocc,occ_c[t,:,nbc],width=0.5,bottom=eigenvalues[:,nbands-1],color='red',edgecolor='none')

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

