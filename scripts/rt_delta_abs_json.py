import matplotlib
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

parser = argparse.ArgumentParser(description='Map of a double-grid')
parser.add_argument('-f' ,'--folder'    , help='Folder with real-time simulations')
parser.add_argument('-j' ,'--job'       , help='Name of job (ex: DELTA-1E+03)')
args = parser.parse_args()

print "Folder: ",args.folder
print "Job: ", args.job

folder = args.folder
job = args.job

prefix = 'B-' # (B)SE, (IP)
print 'Prefix: "B-"'

#pack_files_in_folder(folder)

data = YamboAnalyser(folder)
output = data.get_data(job)

# keys to read the outputs in order
keys=[]
for key in output:
    if key.startswith(prefix+job):
        keys.append(key)

keys=sorted(keys)
print keys

# times to print in file
times=[]
for i in range(0,len(keys)):
    t=re.search("[t]\d{1,}", keys[i]).group()
    print t
    times.append(t)

print times

# output[key][n][0] is x
# output[key][n][1] is y

# Number of lines (number of x points)
nlines = len(output[keys[0]])
# We want to create a file that is x, y1, y2, ..., yn
array = np.zeros((nlines,len(keys)+1))
diff = np.zeros((nlines,len(keys)))

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
        diff[l][i]=y1-y0


# Writing

file1=open(job+'.dat','w')
string = 'eV'
for t in times:
    string = string + '\t' + t
np.savetxt(file1,array,delimiter='\t',header=string)
file1.close()

file2=open(job+'_diff.dat','w')
string = 'eV'
for i,t in enumerate(times):
    if i==0:
       continue
    string = string + '\t' + t + '-t0'
np.savetxt(file2,diff,delimiter='\t',header=string)
file2.close()

print 'Writing done.'

# Plotting
for i in range(1,len(diff[0])):
    plt.plot(diff[:,0],diff[:,i],'-',label='%s-%s'%(times[i],times[0]))
plt.legend(loc=2)
plt.xlabel("eV")
plt.ylabel("TA (arb. units)")
plt.title("TA, job %s in folder %s"%(job,folder))
plt.show()


print 'Done.'
