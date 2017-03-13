import os
from sys import argv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import loadtxt, reshape, zeros, interp, shape
from pylab import *
from matplotlib import colors

# Plot configuration
plt.rc('text', usetex=True)
plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
rcParams['axes.linewidth'] = 2
####################################################################
# Files
prefix   = 'QSSIN-D-100.0fs-1.94eV-300K-DG'
folder   = 'rt-24x24'
path     = folder+'/'+prefix
# use occ_bands.in
occ_file = '%s/%s/o-pulse.YPP-RT_occ_bands_iT' % (folder,prefix)
ext_file = '%s/%s/pulse/o-pulse.external_field'      % (folder,prefix)
name_aux = 'band-occ'
####################################################################
# Variables
t_range        = range(1,101)
ntime          =     101
t_initial      =     1
t_final        =     4
t_step         =     1
t_length       =   5000
t_carrier      =    50    # carrier time
t_output       =     1    #  output time
####################################################################
nbands         =  4
occ_holes      =  1
occ_electrons  =  1
####################################################################
ext = loadtxt(ext_file)
field   = ext[:,2]/max(abs(ext[:,2]))
#fluence = ext[:,8]/max(abs(ext[:,8]))

# Energy range
emin, emax = -1.5, 4.0


##### For energy distribution
##################################
nblines_fit = 1000 # nb of lines for fit files
nblines_data = 8760 # nb of lines for DATA file
#col_nb = 17 # 17 cols per file : #1 is energy, #2-16 is data

#number of files
nb_times = t_length/t_carrier+1
fnb = int(math.ceil(nb_times / 16.0)) # col_nb would need to be there

xval_data = np.zeros([nblines_data]) # only contains the energy values (Data)
xval_fit = np.zeros([nblines_fit,2]) # only contains the energy values (Fit, e and h)
data = np.loadtxt('%s/o-pulse.YPP-RT_occupations_FIT_electrons_1_of_%s'%(path,fnb))
xval_fit[:,0] = data[:,0]
data = np.loadtxt('%s/o-pulse.YPP-RT_occupations_FIT_holes_1_of_%s'%(path,fnb))
xval_fit[:,1] = data[:,0]
data = np.loadtxt('%s/o-pulse.YPP-RT_occupations_DATA_1_of_%s'%(path,fnb))
xval_data[:] = data[:,0]


occ_energy = np.zeros([nblines_data,nb_times]) # ind1: line # (eV value) | ind2: carrier occ at time t_n (both e and h)
fit_energy = np.zeros([nblines_fit,nb_times,2]) # ind1: line # (eV value) | ind2: carrier occ at time t_n | ind3: e/h

# temp arrays to load data
temp_e = np.zeros([16*fnb])
temp_h = np.zeros([16*fnb])
data = np.zeros([nblines_data,17])
fit_e = np.zeros([nblines_fit,17])
fit_h = np.zeros([nblines_fit,17])
#####################################
#####################################

gs = gridspec.GridSpec(8,2)
#fig = plt.figure(figsize=(6,5))
#ax  = fig.add_subplot(111)
#fig.subplots_adjust(left=0.2)

# Plot the band structure
data_lda    = loadtxt('rt-24x24/QSSIN-100.0fs-1.94eV-300K/o-pulse.YPP-RT_occ_bands_iT1')
nkpoints    = shape( data_lda )[0]/nbands
kx          = reshape( data_lda[:,0], (nbands,nkpoints) )[0,:]
knew        = reshape( data_lda[:,0], (nbands,nkpoints) )
elda        = reshape( data_lda[:,1], (nbands,nkpoints) )
nk          = data_lda.shape[0]/nbands
#pos_grid    = [ kx[0], kx[-1]*1./3, kx[-1]/2,kx[-1]*2./3, kx[-1]]
pos_grid    = [ kx[0], kx[-1]*1./4, kx[-1]*5.0/8.6,kx[-1]]
label_grid  = ['$\Gamma$', 'M','K','$\Gamma$']
#label_grid  = ['', 'K','M','K\'','']
#short_grid        = [kx[0], kx[-1]*1./3, kx[-1]/2, kx[-1] ]
#short_label_grid  = ['','K','M','K\'','']

# Plot the occupancies
colors = ['Blue','Blue','Red','Red']
rn     = [ occ_holes, occ_holes, occ_electrons, occ_electrons ]

occ_bands = zeros([ntime,nbands,nk])
for time in t_range:
  for ib in range(nbands):
    abs( reshape( loadtxt( '%s%d' % ( occ_file,time) )[:,2], (nbands,nk) )[ib] )
    occ_bands[time-t_range[0],ib,:] = abs( reshape( loadtxt( '%s%d' % ( occ_file,time) )[:,2], (nbands,nk) )[ib] )

max_occ = max(occ_bands.flatten())

os.system('mkdir -p occupations')
os.system('mkdir -p occupations/%s'%prefix)


####### Fill tables
#fid = 0
#it = 1
##for time in range(nb_times): # could be optimized to only do at time in time_range
#for time in time_range:
#    print time
#    # the first time in each file has time%16=0
#    # eq. to first time/16 = t_carrier
#    #if time%16.0 == 0: # need to change file
#    if (time/16.0)/t_carrier+1 != fid
#        fid += (time/16.0)/t_carrier+1
#        print fid
#        data = np.loadtxt('%s/o-pulse.YPP-RT_occupations_DATA_%s_of_%s'%(path,fid,fnb))
#        fit_e = np.loadtxt('%s/o-pulse.YPP-RT_occupations_FIT_electrons_%s_of_%s'%(path,fid,fnb))
#        fit_h = np.loadtxt('%s/o-pulse.YPP-RT_occupations_FIT_holes_%s_of_%s'%(path,fid,fnb))
#        it = 1
#
#    occ_energy[:,time] = data[:,it]
#    fit_energy[:,time,0] = fit_e[:,it]
#    fit_energy[:,time,1] = fit_h[:,it]
#    it += 1



########
fid = 0
it = 1
for time in t_range:

  # ARRAY FILLING
  print time
  #print "Test value for fid : ",int(time/16.0)+1,(time/16.0)/t_carrier+1
  # the first time in each file has time%16=0
  # eq. to first time/16 = t_carrier
  #if time%16.0 == 0: # need to change file
  if time/16+1 != fid:
      fid =time/16+1
      print fid

      # Temperature data
      f = open('%s/o-pulse.YPP-RT_occupations_FIT_electrons_%s_of_%s'%(path,fid,fnb))
      for line in f:
        if "T(e)" in line:
          lste = map(float,line[10:].strip().split())
          while len(lste)<16:
              lste.append(0)
          temp_e[16*(fid-1):fid*16] = lste
      f.close()
      f = open('%s/o-pulse.YPP-RT_occupations_FIT_holes_%s_of_%s'%(path,fid,fnb))
      for line in f:
        if "T(h)" in line:
          lste = map(float,line[10:].strip().split())
          while len(lste)<16:
              lste.append(0)
          temp_h[16*(fid-1):fid*16] = lste
      f.close()
      #####

      # occupation data
      data = np.loadtxt('%s/o-pulse.YPP-RT_occupations_DATA_%s_of_%s'%(path,fid,fnb))
      fit_e = np.loadtxt('%s/o-pulse.YPP-RT_occupations_FIT_electrons_%s_of_%s'%(path,fid,fnb))
      fit_h = np.loadtxt('%s/o-pulse.YPP-RT_occupations_FIT_holes_%s_of_%s'%(path,fid,fnb))
      it = 1

  occ_energy[:,time] = data[:,it]
  fit_energy[:,time,0] = fit_e[:,it]
  fit_energy[:,time,1] = fit_h[:,it]
  it += 1
  #############

  fig = plt.figure(figsize=(10,6))
  ax  = fig.add_subplot(111)
  fig.subplots_adjust(left=0.2)
  print 'slide number ', time

  # Electron band and occupation window

  b1 = plt.subplot(gs[0:6,0])
  #plt.axis([0.35, kx[-1]*5./6., emin, emax])
  for iband in range(nbands):
    size = occ_bands[time-t_range[0],iband,:]*rn[iband]
    #plt.bar( occ[iband*nk:(iband+1)*nk,0], size, width=0.02, bottom=occ[iband*nk:(iband+1)*nk,1], color=colors[iband], hold=None,linewidth=0,edgecolor='none')
    #plt.bar( kpath[time,iband,:], size, width=0.02, bottom=eband[time,iband,:], color=colors[iband], hold=None,linewidth=0,edgecolor='none')
    plt.bar( knew[iband,:], size, width=0.02, bottom=elda[iband,:], color=colors[iband], hold=None,linewidth=0,edgecolor='none')
    #b1.scatter(occ[iband*nk:(iband+1)*nk,0],occ[iband*nk:(iband+1)*nk,1],s=size,facecolor=colors[iband],edgecolor='none')

  txt = text(0.10,0, '%d fs'%( ext[(time-1)*t_carrier/t_output,0] ),color='k',fontsize=15)
  #txt = text(4,1, '%d fs'%( (time-1)*t_bands ),color='k',fontsize=15)

  b1 = plt.subplot(gs[0:6,0])
  #b1.patch.set_facecolor('MidnightBlue')
#  b1.patch.set_facecolor('white')
#  b1.patch.set_alpha(1.0)
  plt.axis([0.3, kx[-1]*5./6., emin, emax])
  #plt.xticks(short_grid, short_label_grid)
  plt.xticks(pos_grid, label_grid)
  plt.ylabel('E (eV)')
  for band in elda:
    plt.plot(kx, band,'-',color='k',lw=1.0)
  axvline(pos_grid[1] ,color='k')
  axvline(pos_grid[2] ,color='k')
  axvline(pos_grid[3] ,color='k')
  #axvline(kx[-1]*1./6,color='k')
  #axvline(kx[-1]*2./6,color='k')

  #dots.remove()
  name = 'occupations/'+prefix+'/%s%06d' % (name_aux, time)


#  plt.axis([0 ,t_length,-1.2,1.2])
#  plt.xticks([],[])
#  plt.yticks([],[])
#  plt.ylabel('Fluence',fontsize=15)
#  plt.xlim((0,t_length))
#  b2.plot(ext[0:(time-1)*t_bands/t_pulse,0],fluence[0:(time-1)*t_bands/t_pulse],'m-',lw=2)

#  b3 = b2.twinx()

  # External field window
  b2 = plt.subplot(gs[7,0])
  plt.xticks([],[])
  plt.yticks([],[])
  plt.ylabel('Field',fontsize=15)
  plt.xlim((0,t_length))
  plt.ylim((-1.2,1.2))
  plt.plot(ext[0:(time-1)*t_carrier/t_output,0],field[0:(time-1)*t_carrier/t_output],'g-',lw=2)
  if time==t_final:
    plt.plot(ext[:,0],ext[:,2],'g-',lw=2)

  # Plot for the electrons, top right
  b3 = plt.subplot(gs[0:3,1])
  plt.plot(xval_data[:],occ_energy[:,time],'k*',xval_fit[:,0],fit_energy[:,time,0],'r-')
  plt.xlim((1.4,2.5))
  plt.ylim((-0.001,0.02))
  #plt.text(1,1,'T$_e$=%dK'%(temp_e[time]))
  plt.text(0.8,0.9,'T$_e$=%dK'%(temp_e[time]),transform = ax.transAxes)

  # Plot for the holes, bottom right
  b4 = plt.subplot(gs[4:8,1])
  plt.plot(xval_data[:],occ_energy[:,time],'k*',xval_fit[:,1],fit_energy[:,time,1],'b-')
  plt.xlim((-2.2,-1.0))
  plt.ylim((-0.02,0.001))
  plt.text(0.8,0.1,'T$_h$=%dK'%(temp_h[time]),transform = ax.transAxes)


  #plt.show()
  savefig( name ,transparent=False,dpi=300)
  plt.close()
  b3 = b2.twinx()

#plt.show()
# Creating the gif movie
#os.system('convert -delay 40 occupations/%s/%s* %s.gif' % (prefix,name_aux, prefix))
