import os
import imageio
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob, os
from numpy import loadtxt, reshape, zeros, interp, shape
from pylab import *
from matplotlib import colors
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-p", "--path", dest="path",required=True,
                    help="read data from PATH", metavar="PATH")
parser.add_argument("-j", "--job", dest="job",required=True,
                    help="jobname (to be used to select filename)")
parser.add_argument("-f", "--file", dest="file",required=True,
                    help="out filename")
parser.add_argument("-t ", "--time", dest="t_range",nargs='+',help="Time range (units of delta_T)",default=-1,type=int)
parser.add_argument("-g ", "--gif", help="Animate at the end",action='store_true')

args = parser.parse_args()

# Plot configuration
plt.rc('text', usetex=True)
plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
rcParams['axes.linewidth'] = 2
#
occ_file = 'o-'+args.job+'.YPP-RT_occ_bands_iT' 
ntime=0
bands=[]
times=zeros([1000])
for file in sorted(glob.glob(args.path+"/*"+args.job+".YPP-RT_occ_bands_iT*")):
    ntime+=1
    i_time=int(file.split("_iT")[1])
    file_obj = open(file, "r")
    for line in file_obj: 
      if "# TIME=" in line: 
         ch=line.split("# TIME=")[1]
         ch=ch.replace("fs","")
         times[i_time-1]=float(ch)
      if ntime == 1:
        if "|k" in line: 
         nk=0
         label_name = []
         label_pos=[]
         nsymm=0
        if "# Time-window" in line: 
         ch=line[3:].split("#")[0]
         b_ch=ch.split("|")
         t_start=float(b_ch[0])
         t_end=float(b_ch[1])
        if "# Number of bands" in line: 
         ch=line[3:].split("#")[0]
         b_ch=ch.split("|")
         bands.append(int(b_ch[0]))
         bands.append(int(b_ch[1]))
         nbands=int(b_ch[1])-int(b_ch[0])+1
        if "# | TimeStep=" in line: 
         b_ch=line.split("# | TimeStep=")
         t_step=float(b_ch[1].split(" ")[1])
        if not "#" in line: 
         #print line,nk
         nk+=1
         if "[" in line: 
          b_ch=line.split("[")[1]
          nsymm+=1
          ch=b_ch.replace("]","")
          ch=ch.replace("\n","")
          label_name.append(ch)
          label_pos.append(nk)
    file_obj.close() 

#
print "TIME steps    [read]",ntime
print "TIME range[fs][read]",t_start,t_end
if(isinstance(args.t_range,list)):
 print "TIME range[fs][used]",float(args.t_range[0])*t_step+t_start,float(args.t_range[1])*t_step+t_start
print "TIME points   [read]",times[:ntime]
print "K-points           ",nk
print "Bands              ",bands
print "High Symm points   ",label_name
print "          positions",label_pos
#
if(isinstance(args.t_range,list)):
 t_range= args.t_range
 ntime=args.t_range[1]-args.t_range[0]+1
 t_range        = range( args.t_range[0], args.t_range[1]+1,1)
else:
 t_range        = range(1,ntime+1,1)
#
# Energy range
emin, emax = -5.0, 3.0
#
gs = gridspec.GridSpec(8,1)
#fig = plt.figure(figsize=(6,5))
#ax  = fig.add_subplot(111)
#fig.subplots_adjust(left=0.2)

# Plot the occupancies (select the color for each band
#occ_holes      =  1
#occ_electrons  =  1
#colors = ['Blue','Blue','Blue','Blue','Blue','Blue','Red','Red','Blue','Blue','Red','Red']
#rn     = [ occ_holes, occ_holes, occ_electrons, occ_electrons, occ_holes, occ_holes, occ_electrons, occ_electrons, occ_holes, occ_holes, occ_electrons, occ_electrons   ]

occ_bands = zeros([max(t_range),nbands,nk])
ene_bands = zeros([nbands,nk])
max_occ=0
sys.stdout.write("I/O [%s]" % (" " * ntime))
sys.stdout.flush()
sys.stdout.write("\b" * (ntime+1)) # return to start of line, after '['
for file in sorted(glob.glob(args.path+"/*"+args.job+".YPP-RT_occ_bands_iT*")):
    i_time=int(file.split("iT")[1])
    if (i_time>max(t_range) or i_time<min(t_range)) : continue
    file_obj = open(file, "r")
    ib=0
    nk=0
    for line in file_obj: 
     if "|k" in line: 
      nk=0
      ib+=1
     if not "#" in line: 
       nk+=1
       b_ch=line.split()
       occ_bands[i_time-1,ib-1,nk-1] = float(b_ch[2])
       if (abs(occ_bands[i_time-1,ib-1,nk-1])>max_occ): max_occ=occ_bands[i_time-1,ib-1,nk-1]
       ene_bands[ib-1,nk-1] = float(b_ch[1])
    file_obj.close() 
    sys.stdout.write("#")
    sys.stdout.flush()
sys.stdout.write("\n")

if (max_occ == 0): 
 print "Zero occupations. Nothing to plot"
 sys.exit()

images=[]
for i_time in t_range:
  plt.clf()
  sys.stdout.write("Slide#"+`i_time`+"-bands [%s]" % (" " * nbands))
  sys.stdout.flush()
  sys.stdout.write("\b" * (nbands+1)) # return to start of line, after '['
  #
  time=times[i_time-1]
  fig = plt.figure(figsize=(6,5))
  ax  = fig.add_subplot(111)
  fig.subplots_adjust(left=0.2)
  #print 'Slide number @', time,'fs. Band...',

  # Electron band and occupation window

  b1 = plt.subplot(gs[0:6])
  plt.title("Time:  "+str(time)+" fs")
  #plt.axis([0.35, kx[-1]*5./6., emin, emax])
  for iband in range(nbands):
    size = occ_bands[i_time-1,iband-1,:]/max_occ
    if (ene_bands[iband-1,0]<0): 
      band_color='Blue'
    else:
      band_color='Red'
    plt.bar( range(nk), size, width=1, bottom=ene_bands[iband-1,:], color=band_color, hold=None,linewidth=0,edgecolor='none')
    sys.stdout.write("#")
    sys.stdout.flush()
  sys.stdout.write("\n")

  sys.stdout.write("Slide#"+`i_time`+"-plots [%s]" % (" " * nbands))
  sys.stdout.flush()
  sys.stdout.write("\b" * (nbands+1)) # return to start of line, after '['
  b1 = plt.subplot(gs[0:6])
  plt.axis([0, nk-1, emin, emax])
  plt.xticks(label_pos, label_name)
  plt.ylabel('E (eV)')
  for ix in label_pos:
      axvline(ix ,color='k')
  for iband in range(nbands):
      #print 'band (plot)', iband
      plt.plot(range(nk), ene_bands[iband-1,:],'-',color='k',lw=1.0)
      sys.stdout.write("#")
      sys.stdout.flush()
  sys.stdout.write("\n")

  #os.chdir(args.path)
  name = '%s%06d' % (args.path+'/'+args.file, i_time)
  name = './'+name+'.png'
  print 'name ', name

  savefig( name ,transparent=False,dpi=300)
  images.append( imageio.imread(name) )

if (args.gif) :
 filename = args.file + ".gif"
 imageio.mimsave(args.path+'/'+filename, images,fps=1)

