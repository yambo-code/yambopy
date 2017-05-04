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

### aliases
times = [i * 1e15 for i in yrt.times] # carriers output times, in fs
nbands = yrt.nbands # number of bands in the RT simulation

if nbv+nbc != nbands:
    raise NameError('Incompatible number of bands, set nbv and nbc in script.')

## 'path-plot' variables
kindex = yrt.bands_indexes # kpoint indexes (in order) to draw path
eigenvalues = yrt.eigenvalues[kindex,:] # eigenvalues of the bands included in the RT simulation
#
max_occ = np.amax(yrt.occupations[:,kindex,:]) # used to size the distribution plots
occupations = yrt.occupations[:,kindex,:]/max_occ*occ_scaling # format time,kindex,band index (from 0 to nbands, only on path)
norm=Normalize(vmin=0, vmax=occ_scaling, clip=False) # normalizatin class for the color gradiant on bands
#
xocc = np.arange(len(kindex)) # array of ints to plot occupation on path properly
##

## 'fit' variables and function
# FD distrib for fit
def fermi_dirac(E,a,T): # declare E first for fit
    return 1/(1+np.exp((E-a)/T))
#
KtoeV = 8.61733e-5
#
# xeng is an array of values to plot the fit properly
xeng = np.linspace(np.amin(eigenvalues[:,list(range(nbv))]), np.amax(eigenvalues[:,list(range(nbv,nbands))]),1000)
##

##############
# EXT. FIELD #
##############

# The external field is read from the o- file
ext = np.loadtxt('%s/%s/pulse/o-pulse.external_field'%(folder,calc))
field = ext[:,2]/max(abs(ext[:,2])) # polarization : x=1,y=2,z=3

##################
# ENERGY DISTRIB #
##################

# Sort the (n,k) pairs between positive and negative energies
# (If the same energy appears twice, it must not be summed over)
list_e=[] ; list_h=[]
for k in range(yrt.nkpoints):
    for n in range(yrt.nbands):
        e = yrt.eigenvalues[k,n]
        if e<=0.0:
            list_h.append((k,n))
        else:
            list_e.append((k,n))


# Build the occupation tables occ_x[t,(nk)_index,(e|occ)]

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


# *(-1) on holes to fit the same way as electrons
occ_h *= -1



#################
# BAR PLOT DATA #
#################

# occupations in CBs/VBs are summed for easier reading if there are more than one
# Recall that 'occupations' was normalized then multiplied by occ_scaling (for esthetics)
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



####################
# TIME LOOP & PLOT #
####################

# Gridspec allows to place subplots on a grid
# spacing for exemple can be customised
gs = gridspec.GridSpec(9, 8)
#gs.update(hspace=1.1) # bigger horizontal spacing


for t in range(len(times)):
    i=t
    print times[i]

    plt.figure().suptitle('Occupation of the bands and fit to the Fermi-Dirac distribution',fontsize=14,ha='center')

####### bandstructure w/ occupation plot

    ax1 = plt.subplot(gs[0:-1,0:-2])
    # remove x ticks
    ax1.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')

    ax1.set_xlim((0,xocc[-1]))

    # Plot band structure
    ax1.plot(eigenvalues,'k-',lw=2,zorder=0)

    ## Colored lines when degen is beyond degen_thres
    # For that, we compare eigenvalues of (k,n) with (k,n+1)
    # eigenvalues[kindex,n], nbc, nbv
# TODO : loop using parameters
    # VB
    diff_eigen = abs(eigenvalues[:,0]-eigenvalues[:,1])
    ax1.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,0],s=30, c=occ_v[t,diff_eigen>degen_thres,0],cmap='plasma',alpha=1,edgecolors='none',norm=norm)
    ax1.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,1],s=30, c=occ_v[t,diff_eigen>degen_thres,1],cmap='plasma',alpha=1,edgecolors='none',norm=norm)
    # CB
    diff_eigen = abs(eigenvalues[:,2]-eigenvalues[:,3])
    ax1.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,2],s=30, c=occ_v[t,diff_eigen>degen_thres,0],cmap='bwr',alpha=1,edgecolors='none')
    ax1.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,3],s=30, c=occ_v[t,diff_eigen>degen_thres,0],cmap='bwr',alpha=1,edgecolors='none')

    ## occupation in the form of histograms
    # small y-shift for better reading
    ax1.bar(xocc,occ_v[t,:,nbv],width=0.4,bottom=eigenvalues[:,nbv-1]+0.1,color='blue',edgecolor='none')
    ax1.bar(xocc,occ_c[t,:,nbc],width=0.4,bottom=eigenvalues[:,nbands-1]+0.1,color='red',edgecolor='none')

    # text and labels
    ax1.set_ylabel('E (eV)')
    ax1.text(0.10,0.05, '%d fs'%times[t],size=16,transform=ax1.transAxes)

######## field plot

    ax2 = plt.subplot(gs[-1,:])
    ax2.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    ax2.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
    ax2.set_ylabel('Field')
    ax2.set_xlim((0,times[-1]))
    ax2.set_ylim((-1.3,1.3))
    ax2.plot(field[:int(times[t])])

## Rotated version
    ax3 = plt.subplot(gs[0:4,-2:])
    ax4 = plt.subplot(gs[4:8,-2:])
    ax3.text(0.98,0.5,'Energy (eV)',rotation='vertical')
#
## Make the spacing between the two ax3es a bit smaller
#    plt.subplots_adjust(wspace=0.1)
#
# plot the same data on both ax3es
    try: # does not break if fit is not found
        fit,cov = curve_fit(fermi_dirac,occ_e[i,:,0],occ_e[i,:,1])
    except RuntimeError:
        fit=np.array([0,0])

    ax3.scatter(occ_e[i,:,1],occ_e[i,:,0],color='black')
    ax3.plot(fermi_dirac(xeng,fit[0],fit[1]),xeng,'r-')
    ax3.text(0.5,0.9,'Electrons\nT = %d K'%(fit[1]/KtoeV),transform=ax3.transAxes,ha='center',va='center')

    try:
        fit,cov = curve_fit(fermi_dirac,occ_h[i,:,0],occ_h[i,:,1])
    except RuntimeError:
        fit=np.array([0,0])

    ax4.scatter(occ_h[i,:,1],-occ_h[i,:,0],color='black')
    ax4.plot(fermi_dirac(xeng,fit[0],fit[1]),-xeng,'b-')
    ax4.text(0.5,0.1,'Holes\nT = %d K'%(fit[1]/KtoeV),transform=ax4.transAxes,ha='center',va='center')

# zoom-in / limit the view to different portions of the data
    ax3.set_ylim(min(occ_e[i,:,0])-0.1,max(occ_e[i,:,0])+0.1)
    ax4.set_ylim(min(-occ_h[i,:,0])-0.1,max(-occ_h[i,:,0])+0.1)

    ax4.set_xlim(-0.1*max_occ,1.1*max_occ)
    ax3.set_xlim(-0.1*max_occ,1.1*max_occ)

# hide some ticks/labels
    ax3.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    ax3.yaxis.tick_right()

    ax4.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    ax4.yaxis.tick_right()

    plt.show()










