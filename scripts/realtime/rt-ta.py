#import matplotlib
#matplotlib.use('Agg') # prevents crashes if no X server present
from yambopy import *
import sys
import os
import argparse
import re
import numpy as np
import matplotlib.pyplot as plt

"""
After real-time calculations are completed, plots the Transient Absorption (TA) at different times.
Setting inside the script : prefix
"""

parser = argparse.ArgumentParser(description='Make Transiant Absorption (TA) plots from RT-BSE')
parser.add_argument('-f' ,'--folder'    , help='Real-time folder (e.g.: rt-24x24)')
parser.add_argument('-j' ,'--job'       , help='Name of job (e.g.: QSSIN-2.0eV)')
parser.add_argument('-p' ,'--prefix'    , help='Prefix of the BSE calculations (e.g.: B-XRK-XG)')
parser.add_argument('-nt','--notext'    , help='Skips the writing of the data', action='store_false')
parser.add_argument('-np','--nopack'    , help='Skips packing o- files into .json files', action='store_false')
args = parser.parse_args()

folder = args.folder
job = args.job
prefix = args.prefix

source = folder+'/'+job

if args.nopack: # True by default, False if -np used
    print 'Packing relevant calculations'
    pack_files_in_folder(source,mask=prefix)
    print 'Done.'

data = YamboAnalyser(source)
output = data.get_data((prefix,'eps'))

# keys to read the outputs in order
keys=output.keys()

# sorting
s=[]
for i,key in enumerate(keys):
        s.append((int(re.search("\d{1,}", key).group()),key))

s=sorted(s) # s is list of tuples (int,str)
            # having integers allows to sort correctly

times = ['t'+str(i[0]) for i in s]
keys  = [i[1] for i in s]

print "Sorted keys: ",keys

# output[key][n][0] is x
# output[key][n][1] is y

# Number of lines (number of x points)
nlines = len(output[keys[0]])
# We want to create a file that is x, y1, y2, ..., yn
array = np.zeros((nlines,len(keys)+1))
diff = np.zeros((nlines,len(keys)))
ymax=0;ymin=0

for l in range(0,nlines):
    # first, value of x
    x = output[keys[0]][l][0]
    array[l][0]=x
    diff[l][0]=x
    # then all the y's
    # t0
    y0 = output[keys[0]][l][1]
    array[l][1] = y0
    # additional columns
    for i,key in enumerate(keys):
        # t_i > t0
        y1 = output[key][l][1]
        array[l][i+1]=y1

        if i==0:
            continue
        diff[l][i]=(y1-y0)/y0
        if max(diff[l][1:])>ymax:
            ymax=max(diff[l][1:])
        if min(diff[l][1:])<ymin:
            ymin=min(diff[l][1:])

ymin=1.1*ymin ; ymax=1.1*ymax


# Writing
if args.notext:
    print 'Writing data to files...'
    f=open(job+prefix+'.dat','w')
    string = 'eV'
    for t in times:
        string = string + '\t' + t
    np.savetxt(f,array,delimiter='\t',header=string)
    f.close()

    f=open(job+'_diff.dat','w')
    string = 'eV'
    for i,t in enumerate(times):
        if i==0:
           continue
        string = string + '\t(' + t + '-t0)/t0'
    np.savetxt(f,diff,delimiter='\t',header=string)
    f.close()

    print 'Writing done.'
else:
    print '"--notext" flag'

# Plotting
os.system('mkdir -p TA/%s'%job)
fig,ax1=plt.subplots()
for i in range(1,len(diff[0])):
    plt.plot(diff[:,0],diff[:,i],'-',label='(%s-%s)/%s'%(times[i],times[0],times[0]))
    plt.fill_between(diff[:,0],diff[:,i],where=diff[:,i]>0,color='r',alpha=0.4)
    plt.fill_between(diff[:,0],diff[:,i],where=diff[:,i]<0,color='b',alpha=0.4)
    plt.legend(loc=2)
    plt.xlabel("eV")
    plt.ylabel("TA (arb. units)")
    plt.axvline(1.94,color='k')
    plt.title("TA, job %s/%s in folder %s"%(job,prefix,folder))
    plt.ylim((ymin,ymax))



    # inset with BSE
    a = plt.axes([0.2, .15, .25, .25])
    plt.plot(array[:,0],array[:,i],'k-',label='BSE@t0')
    plt.axvline(1.94,color='k')
    plt.savefig('TA/%s/%d'%(job,i))
    plt.close()


#os.system('convert -delay 40 TA/%s/*.png TA_%s.gif'%(job,job) )


print 'Done.'
