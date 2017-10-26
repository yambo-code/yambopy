from __future__ import print_function
# Version : dec 2nd
from builtins import str
from builtins import range
import os
import re
import numpy as np

# Get a list of all exciton files, make a new file with x=time, y1 = exc A1, y2 = exc A2

# List of exciton files
dirfiles = os.listdir('.')
files=[]
for f in dirfiles:
    if '_E_' in f:
        files.append(f)

files.sort()
print(files)

# Get times
times=[]
for f in files:
    time=re.search('[0-9]{3,}',f).group()
    times.append(time)

# Get excitons
excitons=[]
i_thr=0.005 # intensity threshold
exc_list = (1,2) # excitons wanted (here 1st and 2nd one that respect the threshold)

for p,f in enumerate(files):
    n=1
    print(files[p],times[p])
    content = np.loadtxt(f)
    for e,i,index in content:
        if i>i_thr:
            if n in exc_list:
                n=n+1
                excitons.append((times[p],e))
# excitons is an array of the form
# excitons[(time0,energy_n),(time0,energy_n+1),...,(timep,energy_n+m)]


f = open('output.dat','w')
for i in range(0,len(excitons)):
    t = excitons[i][0]
    e = excitons[i][1]
    line = str(t)+'\t'+str(e)+'\n'
    f.write(line)

f.close()
