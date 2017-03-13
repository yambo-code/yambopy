import os
from matplotlib.cm import ScalarMappable
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle, Polygon, Circle, Arrow, FancyArrowPatch
from matplotlib import colors, colorbar
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import loadtxt, reshape, zeros, interp, shape, arange, max, sqrt, min
from pylab import *
from matplotlib import colors
import matplotlib.patches as patches
import itertools
from yambopy import *

# Plot configuration
plt.rc('text', usetex=True)
plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
rcParams['axes.linewidth'] = 2
####################################################################
# Files
#prefix   = 'QSSIN-D-dostest-1.94eV-300K-0.05fs-DG'
prefix   = 'RTstep-2.0'
#folder   = 'rt-24x24/QSSIN-D-dostest-1.94eV-300K-0.05fs-DG'
folder   = 'rt-24x24/RTstep-2.0'
# use occ-t.in
occ_file = '%s/o-pulse.YPP-RT_occupations_'    % (folder)
fluence  = loadtxt('%s/pulse/o-pulse.external_field' % (folder) )
#kerr     = loadtxt('kerr-24x24/ip-%s.dat'   % (prefix) )
####################################################################
nbands         =  2
nkpoints       =  300
nkfull         =  576
nt             =  11
start_band     =  25
radius_k       = 1000
###################################################################
time_limit  = 500
output_time = 50
#rn = 500 not used
time_string = range(1,11,1)# + range(20,101,10)
#####################################################################
B1 = array( [ 1.026342, 0.592559,  0.000000])
B2 = array( [ 0.000000, 1.185117,  0.000000])
B3 = array( [ 0.000000, 0.000000,  0.157080])
x1 =  1./3.*B1[:2] + 1./3.*B2[:2]
x2 = -1./3.*B1[:2] + 2./3.*B2[:2]
x3 = -2./3.*B1[:2] + 1./3.*B2[:2]
x4, x5, x6 = -x1, -x2, -x3
x7, x8, x9  = x1 + B1[:2], x2 + B1[:2], x3 + B1[:2]
x10,x11,x11 = x4 + B1[:2], x5 + B1[:2], x6 + B1[:2]

hexagon = [ [x1,x2,x3,x4,x5,x6] ]
for ik in [-1,1]:
  hexagon.append( [x1+ik*B1[:2],x2+ik*B1[:2],x3+ik*B1[:2],x4+ik*B1[:2],x5+ik*B1[:2],x6+ik*B1[:2]] )
  hexagon.append( [x1+ik*B2[:2],x2+ik*B2[:2],x3+ik*B2[:2],x4+ik*B2[:2],x5+ik*B2[:2],x6+ik*B2[:2]] )

""" Reading of the data files
"""
print('Reading of the data files...')
k_ibz   = zeros([nkpoints,3])
k_fbz   = zeros([nkfull  ,3])
kpoint  = zeros([nkfull*9,3])
ktim    = zeros([nkfull,nt])
kocc_vb = zeros([nkfull,nbands,nt])
kocc_cb = zeros([nkfull,nbands,nt])

lattice = YamboSaveDB()
count = 0
sym = np.array((-1.0,1.0,1.0))
print len(k_fbz)



for ik in range(1,nkpoints+1):
  #print occ_file + 'k%d_b%d' % (ik,start_band)
  print ik
  ktim[ik-1] = loadtxt( occ_file + 'k%d_b%d' % (ik,start_band) )[:nt,0]
  for ib in range(nbands):
    #print ik,ib
    kocc_vb[ik-1,ib] = abs(loadtxt( occ_file + 'k%d_b%d' % (ik,ib+start_band  ) )[:nt,1])
    kocc_cb[ik-1,ib] = abs(loadtxt( occ_file + 'k%d_b%d' % (ik,ib+start_band+2) )[:nt,1])
    kaux = open(occ_file + 'k%d_b%d' % (ik,ib+start_band), 'r').readlines()
    for line in kaux:
      if '( cc)' in line:
        kref = map(float, line.replace(':','').split()[3:6])
        k_ibz[ik-1] = kref

k_fbz[:nkpoints] = k_ibz
i_fbz = nkpoints
for i_ibz,k_0 in enumerate(k_ibz):
  for op in lattice.sym_car:
    k_proj = np.dot(op,k_0)
    i_1 = 0
    exist = False
    while i_1 < nkfull:
        error_k = abs(np.dot(k_proj - k_fbz[i_1],k_proj - k_fbz[i_1]))
        if error_k < 0.001: # this value already exists, get out of the loop
            exist = True
            break
        i_1 += 1

    if not exist: # if it does not exist, we need to add it
        k_fbz[i_fbz] = k_proj
        kocc_vb[i_fbz] = kocc_vb[i_ibz]
        kocc_cb[i_fbz] = kocc_cb[i_ibz]
        i_fbz += 1
        print i_fbz

for ik4,kref in enumerate(k_fbz):
  ik3 = 0
  for ik1,ik2 in itertools.product(range(3),range(3)):
    kpoint[ nkfull*ik3 + ik4 ] = kref + ik1*B1 + ik2*B2 - (B1+B2)
    ik3 += 1

print('Done!')
max_occ_vb = max(kocc_vb[:,0,:]) + max(kocc_vb[:,1,:])
max_occ_cb = max(kocc_cb[:,0,:]) + max(kocc_cb[:,1,:])
print(max_occ_cb)
print(max_occ_vb)
norm = mpl.colors.Normalize(vmin=0.0, vmax=0.4)
color_map = plt.get_cmap("viridis")

for it in time_string:
  print 'time step', it
  # Occupations and color map
  # Generation of occupation in the replicas of the FBZ
  fbz_vb, fbz_cb = concatenate(9*[kocc_vb[:,0,it]+kocc_vb[:,1,it]]), concatenate(9*[kocc_cb[:,0,it]+kocc_cb[:,1,it]])
  smap = ScalarMappable(cmap=color_map,norm=norm)
  c_vb, c_cb     = smap.to_rgba(sqrt(fbz_vb/max_occ_vb)), smap.to_rgba(sqrt(fbz_cb/max_occ_cb))
  c_vb[:,3], c_cb[:,3] = 1.0, 1.0

  fig    = plt.figure(figsize=(8,8))

  # VB dynamics
  axg_vb = fig.add_axes( [ 0.1, 0.4, 0.4, 0.4 ])
  plt.axis( [-0.8,0.8,-0.8,0.8]  )
  plt.xticks([])
  plt.yticks([])
  #for vertices in hexagon:
  #  axg_vb.add_patch(Polygon(vertices,closed=True,fill=False,color='w',lw=2.5))
  plot_vb=scatter(kpoint[:,0],kpoint[:,1],c=c_vb,marker='H',s=68,edgecolor='none',rasterized=True)
  txt = text(-0.75,0.65, 'VB' ,color='w',fontsize=25)

  # CB dynamics
  axg_cb = fig.add_axes( [ 0.55, 0.4, 0.4, 0.4 ])
  plt.axis( [-0.8,0.8,-0.8,0.8]  )
  plt.xticks([])
  plt.yticks([])
  #for vertices in hexagon:
  #  axg_cb.add_patch(Polygon(vertices,closed=True,fill=False,color='w',lw=2.5))
  plot_cb=scatter(kpoint[:,0],kpoint[:,1],c=c_cb,marker='H',s=68,edgecolor='none',rasterized=True)
  txt2 = text(-0.75,0.65, 'CB' ,color='w',fontsize=25)

  # Field intensity
  axg_kerr = fig.add_axes( [ 0.10, 0.2, 0.85, 0.15 ])
  plt.xticks([])
  plt.yticks([])
  plt.xlim((0,time_limit))
  plt.ylim((0,max(fluence[:,7])*1.05))
  plt.plot(fluence[:it*output_time,0],fluence[:it*output_time,7],'b-',lw=4)
  fill_between(fluence[:it*output_time,0],fluence[:it*output_time,7],0,color='b',alpha=0.5)
  txt3 = text(1000,500, '%d fs'%(ktim[0][it]) ,color='k',fontsize=30)
  #txt4 = text(10, max(fluence[:,7])*0.8, 'Intensity' ,color='b',fontsize=15)
  #txt5 = text(1000,max(fluence[:,7])*0.5, 'Kerr angle' ,color='r',fontsize=15)
  # kerr angle
  k_plot = twinx()
  plt.xticks([])
  plt.yticks([])
  plt.xlim((0,time_limit))
#  plt.ylim((0,max(kerr[:,1])*1.05))
#  plt.plot(kerr[:it,0],kerr[:it,1],'r--',lw=4)
#  fill_between(kerr[:it,0],kerr[:it,1],0,color='r',alpha=0.5)
  #txt5 = text(1300,100, 'Kerr angle' ,color='r',fontsize=20)
  plt.show()
  #plt.show()
  print 'time step'
  print 'movie%06d.png'%it
  savefig('movie%06d.png'%it, transparent=False, dpi=300,bbox_inches='tight')
  txt.remove()
  txt2.remove()
  #txt3.remove()
  #txt4.remove()
  #txt5.remove()
  plot_cb.remove()
  plot_vb.remove()
  plt.close()

os.system('mkdir -p movie/%s'%prefix)
os.system('mv movie*.png movie/%s/'%prefix)
