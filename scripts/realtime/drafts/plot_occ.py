# This script reads the carrier database
# and display it along a path in histogram form

from yambopy import *
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit

############
# SETTINGS #
############

#save = 'rt-24x24/SAVE' # Save with the appropriate band structure
folder = 'rt-24x24'
calc   = 'QSSIN-D-100.0fs-2.07eV-300K-DG' # Where RT carrier output is
path = [[0.0,0.0,0.0],[0.5,0.0,0.0],[0.33333,0.33333,0.0],[0.0,0.0,0.0]]
nbv = 2 ; nbc = 2 # nb of valence and conduction bands
occ_scaling = 1 # max occupation will be 1eV high
degen_thres = 0.1 # Energy below which two bands are considered degenerate


########
# INIT #
########

# Instance containing bandstructure (as used in RT sim) and occupations
yrt = YamboRTDB(folder=folder,calc=calc)

yrt.get_path(path) # Generates kindex

# aliases
kindex = yrt.bands_indexes # kpoint indexes (in order) to draw path
eigenvalues = yrt.eigenvalues[kindex,:] # eigenvalues of the bands included in the RT simulation
max_occ = np.amax(yrt.occupations[:,kindex,:]) # used to size the distribution plots
nbands = yrt.nbands # number of bands in the RT simulation
times = [i * 1e15 for i in yrt.times] # carriers output times, in fs
occupations = yrt.occupations[:,kindex,:]/max_occ*occ_scaling # format time,kindex,band index (from 0 to nbands, only on path)

if nbv+nbc != nbands:
    raise NameError('Incompatible number of bands, set nbv and nbc in script.')

##################
# ENERGY DISTRIB #
##################
#
## Sort the (n,k) pairs between positive and negative energies
#i=0;j=0
#list_e=[] ; list_h=[]
#for k in range(yrt.nkpoints):
#    for n in range(yrt.nbands):
#        e = yrt.eigenvalues[k,n]
#        if e<=0.0:
#            list_h.append((k,n))
#        else:
#            list_e.append((k,n))
#
#
## Build the occupation tables occ_x[t,(nk)_index,(e|occ)]
#
#occ_e = np.zeros((len(times),len(list_e),2))
#for t in range(len(times)):
#    for i,(k,n) in enumerate(list_e):
#        occ_e[t,i,0]=yrt.eigenvalues[k,n]
#        occ_e[t,i,1]=yrt.occupations[t,k,n]
#
#occ_h = np.zeros((len(times),len(list_h),2))
#for t in range(len(times)):
#    for i,(k,n) in enumerate(list_h):
#        occ_h[t,i,0]=yrt.eigenvalues[k,n]
#        occ_h[t,i,1]=yrt.occupations[t,k,n]
#
#
## *(-1) on holes to fit the same way as electrons
#occ_h *= -1
#
#def fermi_dirac(E,a,T): # declare E first for fit
#    return 1/(1+np.exp((E-a)/T))
#
#eVtoK=11599
## TODO : check if degen is summed or not (looks like it when comparing to YPP)
## TODO : print the error on the fit somewhere ? (especially bad at early times)
#
#for i,t in enumerate(times):
#    print times[i]
#    # temperature is in fit[0]*eVtoK
#
#
### Rotated version
#    f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
#    f.suptitle('Occupation of the bands and fit to the Fermi-Dirac distribution')
#    ax.set_ylabel('Electrons')
#    ax2.set_ylabel('Holes')
#    f.text(0.98,0.5,'Energy (eV)',rotation='vertical')
#
## Make the spacing between the two axes a bit smaller
#    plt.subplots_adjust(wspace=0.1)
#
## plot the same data on both axes
#    try: # does not break if fit is not found
#        fit,cov = curve_fit(fermi_dirac,occ_e[i,:,0],occ_e[i,:,1])
#    except RuntimeError:
#        fit=np.array([0,0])
#
#    ax.scatter(occ_e[i,:,1],occ_e[i,:,0],color='black')
#    ax.plot(fermi_dirac(occ_e[i,:,0],fit[0],fit[1]),occ_e[i,:,0],'r+')
#    ax.text(0.75,0.8,'Fitted T: %d K'%(fit[1]*eVtoK),transform=ax.transAxes)
#
#    try:
#        fit,cov = curve_fit(fermi_dirac,occ_h[i,:,0],occ_h[i,:,1])
#    except RuntimeError:
#        fit=np.array([0,0])
#
#    ax2.scatter(occ_h[i,:,1],-occ_h[i,:,0],color='black')
#    ax2.plot(fermi_dirac(occ_h[i,:,0],fit[0],fit[1]),-occ_h[i,:,0],'g+')
#    ax2.text(0.75,0.2,'Fitted T: %d K'%(fit[1]*eVtoK),transform=ax2.transAxes)
#
## zoom-in / limit the view to different portions of the data
#    ax.set_ylim(min(occ_e[i,:,0])-0.1,max(occ_e[i,:,0])+0.1)
#    ax2.set_ylim(min(-occ_h[i,:,0])-0.1,max(-occ_h[i,:,0])+0.1)
#
#    ax2.set_xlim(-0.1*max_occ,1.1*max_occ)
#    ax.set_xlim(-0.1*max_occ,1.1*max_occ)
#
## hide the spines between ax and ax2, hide some ticks/labels
#    ax.yaxis.tick_right()
#    ax.spines['bottom'].set_visible(False)
#    ax.xaxis.tick_top()
#    ax.tick_params(labeltop='off')  # don't put tick labels at the top
#
#    ax2.yaxis.tick_right()
#    ax2.xaxis.tick_top()
#    ax2.tick_params(labeltop='off') # don't put tick labels at the top
#    ax2.spines['top'].set_visible(False)
#
#    plt.show()
#
#exit()


#################
# BAR PLOT DATA #
#################

# occupations in CBs/VBs are summed for easier reading if there are more than one
# Recall that 'occupations' was normalized then multiplied by occ_scaling (for esthetics)
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
        np.add(occ_c[:,:,n],occ_c[:,:,nbc],occ_c[:,:,nbc]) # each time we add the occ of the current band to the total


##############
# EXT. FIELD #
##############

# The external field is read from the o- file
ext = np.loadtxt('%s/%s/pulse/o-pulse.external_field'%(folder,calc))
field = ext[:,2]/max(abs(ext[:,2])) # polarization : x=1,y=2,z=3




# Gridspec allows to place subplots on a grid
# spacing for exemple can be customised
gs = gridspec.GridSpec(8, 1)
gs.update(hspace=1.1) # bigger horizontal spacing

# To plot occ properly, we need an order array of integers
xocc = np.arange(len(kindex))
# the occupation is renormalized between 0 and 1
# then multiplied by occ_scaling for esthetics, so vmax=occ_scaling
norm=Normalize(vmin=0, vmax=occ_scaling, clip=False)
#TODO : Use degen_thres and replot the BS using occ(k)/max_occ for color function
for t in range(len(times)):
#for t in (10,):

    ax1 = plt.subplot(gs[:-1, :]) # bandstructure w/ occupation plot
    plt.xticks(()) # remove x ticks
# Band structure
    plt.plot(eigenvalues,'k-',lw=2,zorder=0)
## Colored lines when degen is beyond degen_thres
# For that, compare eigenvalues of (k,n) with (k,n+1)
# eigenvalues[kindex,n], nbc, nbv
    diff_eigen = abs(eigenvalues[:,0]-eigenvalues[:,1])
    plt.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,0],s=30, c=occ_v[t,diff_eigen>degen_thres,0],cmap='plasma',alpha=1,edgecolors='none',norm=norm)
    plt.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,1],s=30, c=occ_v[t,diff_eigen>degen_thres,1],cmap='plasma',alpha=1,edgecolors='none',norm=norm)
    diff_eigen = abs(eigenvalues[:,2]-eigenvalues[:,3])
    plt.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,2],s=30, c=occ_v[t,diff_eigen>degen_thres,0],cmap='bwr',alpha=1,edgecolors='none')
# occupation in the form of histogram
    plt.bar(xocc,occ_v[t,:,nbv],width=0.2,bottom=eigenvalues[:,nbv-1]+0.1,color='blue',edgecolor='none')
    plt.bar(xocc,occ_c[t,:,nbc],width=0.2,bottom=eigenvalues[:,nbands-1]+0.1,color='red',edgecolor='none')

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




##### Old snippets
#### Distribution plot not rotated
#    fig,(ax,ax2) = plt.subplots(1, 2, sharey=True)
#
## 2 axes to make the broken axis
#    fit,cov = curve_fit(fermi_dirac,occ_e[i,:,0],occ_e[i,:,1])
#    ax2.plot(occ_e[i,:,0],fermi_dirac(occ_e[i,:,0],fit[0],fit[1]),'r+')
#    ax2.scatter(occ_e[i,:,0],occ_e[i,:,1],color='black')
#
#    fit,cov = curve_fit(fermi_dirac,occ_h[i,:,0],occ_h[i,:,1])
#    ax.plot(-occ_h[i,:,0],fermi_dirac(occ_h[i,:,0],fit[0],fit[1]),'g+')
#    ax.scatter(-occ_h[i,:,0],occ_h[i,:,1],color='black')
#
## zoom-in / limit the view to different portions of the data
#    ax2.set_xlim(min(occ_e[i,:,0])-0.1,max(occ_e[i,:,0])+0.1)
#    ax.set_xlim(min(-occ_h[i,:,0])-0.1,max(-occ_h[i,:,0])+0.1)
#
#    ax.set_ylim(-0.1*max_occ,1.1*max_occ)
#    ax2.set_ylim(-0.1*max_occ,1.1*max_occ)
#
## hide the spines between ax and ax2
#    ax.spines['right'].set_visible(False)
#    ax2.spines['left'].set_visible(False)
#    ax.yaxis.tick_left()
#    ax.tick_params(labeltop='off') # don't put tick labels at the top
#    ax2.yaxis.tick_right()
#
## More ticks esthetics
#    ax2.xaxis.tick_top()
#    ax2.yaxis.tick_right()
#    ax2.tick_params(labelright='off')
#
#    ax.xaxis.tick_top()
#
## Make the spacing between the two axes a bit smaller
#    plt.subplots_adjust(wspace=0.15)
#
#    plt.show()
###########

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
