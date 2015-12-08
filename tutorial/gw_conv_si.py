#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using Yambo
#
from __future__ import print_function
from yambopy.inputfile import *
from yambopy.outputfile import *
from yambopy.analyse import *
from pwpy.inputfile import *
from pwpy.outputxml import *
import subprocess
import spur

yambo = 'yambo'
shell = spur.LocalShell()

if not os.path.isdir('database'):
    os.mkdir('database')

#check if the nscf cycle is present
if os.path.isdir('nscf/si.save'):
    print('nscf calculation found!')
else:
    print('nscf calculation not found!')
    exit() 

#check if the SAVE folder is present
if not os.path.isdir('database/SAVE'):
    print('preparing yambo database')
    log_p2y   = shell.run(['p2y'],  cwd="nscf/si.save")
    print(log_p2y.output,file=open("p2y.log","w"))
    log_yambo = shell.run(['yambo'],cwd="nscf/si.save")
    print(log_yambo.output,file=open("yambo.log","w"))
    shell.run('mv nscf/si.save/SAVE database'.split())

if not os.path.isdir('gw_conv'):
    os.mkdir('gw_conv')
    shell.run('cp -r database/SAVE gw_conv'.split())

#create the yambo input file
y = YamboIn('%s -d -g n -V all'%yambo,folder='gw_conv')
QPKrange,_ = y['QPkrange']
y['QPkrange'] = [QPKrange[:2]+[6,10],'']
conv = { 'FFTGvecs': [[10,15,20],'Ry'],
         'NGsBlkXd': [[1,2,5], 'Ry'],
         'BndsRnXd': [[1,10],[1,20],[1,30]] }
y.arguments.append('WFbuffIO')

def run(filename):
    """ Function to be called by the optimize function """
    folder = filename.split('.')[0]
    print(filename, folder)
    shell.run(('yambo -F %s -J %s -C %s 2> %s.log'%(filename,folder,folder,folder)).split(),cwd='gw_conv')

y.optimize(conv,run=run)

#pack the files in .json files
for folder in subprocess.check_output('ls gw_conv/*/o-*',shell=True).splitlines():
    folder = '/'.join(folder.split('/')[:-1])
    y = YamboOut(folder)
    if not y.locked():
        y.pack()
        y.put_lock()

#plot the results using yambmo analyser
y = YamboAnalyser('gw_conv')
print(y)
y.plot_gw('qp')
print('done!')
