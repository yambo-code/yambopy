# This script reads the carrier database
# and display it along a path in histogram form
# along with a representation of the carriers in energy space

from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
from yambopy import *
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit
import os

############
# SETTINGS #
############

folder = 'rt-30x30'
calc   = 'QSSIN-100.0fs-2.08eV-300K-DG' # Where RT carrier output is
path = [[0.0,0.0,0.0],[0.5,0.0,0.0],[0.33333,0.33333,0.0],[0.0,0.0,0.0]]
nbv = 2 ; nbc = 2 # nb of valence and conduction bands
occ_scaling = 1 # max occupation will be 1eV high
degen_thres = 0.1 # Energy below which two bands are considered degenerate

########
# INIT #
########

# For saving pictures
os.system('mkdir -p occupations/%s/%s'%(folder,calc))

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
    return old_div(1,(1+np.exp(old_div((E-a),T))))
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
field = old_div(ext[:,2],max(abs(ext[:,2]))) # polarization : x=1,y=2,z=3

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

# y range for band structure & energy plots
ymin_v= np.amin(eigenvalues[:,:nbv])-0.1
ymin_c= np.amin(eigenvalues[:,nbv:])-0.1

ymax_v= max(np.amax(eigenvalues[:,:nbv])+np.amax(occ_c[:,:,nbv:])+0.1, np.amax(eigenvalues[:,:nbv])+0.1)
ymax_c= max(np.amin(eigenvalues[:,nbv:])+np.amax(occ_c[:,:,nbv:])+0.1, np.amax(eigenvalues[:,nbv:])+0.1)
###

for t in range(len(times)):
    i=t
    print(times[i])

    name = 'occupations/'+folder+'/'+calc+'/%d.png' % (times[t])

    fig = plt.figure()

    fig.suptitle('Occupation of the bands and fit to the Fermi-Dirac distribution',fontsize=14,ha='center')

####### bandstructure w/ occupation plot

    ax1c = plt.subplot(gs[0:4,0:-2])
    ax1v = plt.subplot(gs[4:8,0:-2])
    # remove x ticks
    ax1c.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    ax1v.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')

    # set x range
    ax1c.set_xlim((0,xocc[-1]))
    ax1v.set_xlim((0,xocc[-1]))
    # y range is defined with ax3 and ax4 (they share y axis with ax1)


    # Plot band structure
    ax1v.plot(eigenvalues[:,:nbv],'k-',lw=2,zorder=0)
    ax1c.plot(eigenvalues[:,nbv:],'k-',lw=2,zorder=0)

    ## Colored spots when degen is beyond degen_thres
    # For that, we compare eigenvalues of (k,n) with (k,n+1)
# note : if more than 2 VB/CB, this scatter scheme might not be optimal (e.g. with 1 + 2 degen bands)

    # VB
    for n in range(nbv-1): # we compare n and n+1 <= nbv
        # bool array with condition on degeneracy
        diff_eigen = abs(eigenvalues[:,n]-eigenvalues[:,n+1])
        # plot for points that respect the condition
        ax1v.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,n]  ,s=30, c=occ_v[t,diff_eigen>degen_thres,n]  ,cmap='plasma',alpha=1,edgecolors='none',norm=norm)
        ax1v.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,n+1],s=30, c=occ_v[t,diff_eigen>degen_thres,n+1],cmap='plasma',alpha=1,edgecolors='none',norm=norm)

    # CB
    for n in range(nbc-1):
        diff_eigen = abs(eigenvalues[:,nbv+n]-eigenvalues[:,nbv+n+1])
        ax1c.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,nbv+n]  ,s=30, c=occ_c[t,diff_eigen>degen_thres,n]  ,cmap='plasma',alpha=1,edgecolors='none',norm=norm)
        ax1c.scatter(xocc[diff_eigen>degen_thres],eigenvalues[diff_eigen>degen_thres,nbv+n+1],s=30, c=occ_c[t,diff_eigen>degen_thres,n+1],cmap='plasma',alpha=1,edgecolors='none',norm=norm)

    ## occupation in the form of histograms
    # small y-shift for better reading
    ax1v.bar(xocc,occ_v[t,:,nbv],width=0.4,bottom=eigenvalues[:,nbv-1]+0.1,color='blue',edgecolor='none')
    ax1c.bar(xocc,occ_c[t,:,nbc],width=0.4,bottom=eigenvalues[:,nbands-1]+0.1,color='red',edgecolor='none')

    # text and labels
    fig.text(0.05,0.6,'Energy (eV)',size=16,rotation='vertical')
    fig.text(0.50,0.91, '%d fs'%times[t],size=16)

######## field plot

    ax2 = plt.subplot(gs[-1,:])

    # remove ticks and labels
    ax2.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    ax2.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')

    # text
    ax2.set_ylabel('Field')

    # frame size
    ax2.set_xlim((0,times[-1]))
    ax2.set_ylim((-1.3,1.3))

    ax2.plot(field[:int(times[t])])

## Plot of the occupation as a function of energy (rotated to match the band structure)
    ax3 = plt.subplot(gs[0:4,-2:],sharey=ax1c)
    ax4 = plt.subplot(gs[4:8,-2:],sharey=ax1v)

# plot the data
    try: # does not break if fit is not found
        fit,cov = curve_fit(fermi_dirac,occ_e[i,:,0],occ_e[i,:,1])
    except RuntimeError:
        fit=np.array([0,0])

    ax3.scatter(occ_e[i,:,1],occ_e[i,:,0],s=10,color='black')
    ax3.plot(fermi_dirac(xeng,fit[0],fit[1]),xeng,'r-')
    ax3.text(0.5,0.9,'Electrons\nT = %d K'%(old_div(fit[1],KtoeV)),transform=ax3.transAxes,ha='center',va='center')

    try:
        fit,cov = curve_fit(fermi_dirac,occ_h[i,:,0],occ_h[i,:,1])
    except RuntimeError:
        fit=np.array([0,0])

    ax4.scatter(occ_h[i,:,1],-occ_h[i,:,0],color='black')
    ax4.plot(fermi_dirac(xeng,fit[0],fit[1]),-xeng,'b-')
    ax4.text(0.5,0.1,'Holes\nT = %d K'%(old_div(fit[1],KtoeV)),transform=ax4.transAxes,ha='center',va='center')

# set x and y range
    ax4.set_xlim(-0.1*max_occ,1.1*max_occ)
    ax3.set_xlim(-0.1*max_occ,1.1*max_occ)

    ax3.set_ylim(( ymin_c,ymax_c ))
    ax4.set_ylim(( ymin_v,ymax_v ))

# hide some ticks/labels
    ax3.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    ax3.tick_params(axis='y',labelleft='off',labelright='off')

    ax4.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
    ax4.tick_params(axis='y',labelleft='off',labelright='off')

    plt.savefig( name ,transparent=False,dpi=300)
    print(name)
    #plt.show()
    plt.close(fig)










