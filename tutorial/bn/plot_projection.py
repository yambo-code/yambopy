# Author: Alejandro Molina-Sanchez, Henrique Pereira Coutada Miranda
# Projection of the band states into orbial & atom component 
#
from qepy import *
import argparse
from schedulerpy import oarsub

folder = 'nscf'
prefix = 'bn'
mod_espresso= 'espresso/5.1'
qe_nodes=1
qe_cores=1
walltime="1:0:0"

path = Path([ [[0.0, 0.0, 0.0],'G'],
              [[0.5, 0.0, 0.0],'M'],
              [[1./3,1./3,0.0],'K'],
              [[0.0, 0.0, 0.0],'G']], [50,25,50])

#parse options
parser = argparse.ArgumentParser(description='Test the yambopy script.')
parser.add_argument('-c' ,'--calc',    action="store_true", help='Project orbitals')
parser.add_argument('-r' ,'--run',     action="store_true", help='Project orbitals')
parser.add_argument('-a' ,'--analyse', action="store_true", help='Analyse data')
parser.add_argument('-p' ,'--plot',    action="store_true", help='Analyse data')
args = parser.parse_args()

if args.calc:
    f = open('proj.in','w')
    projwfc = ProjwfcIn(prefix)
    projwfc.write(folder=folder)
    projwfc.run(folder=folder)


if args.run:
    qe = oarsub(nodes=qe_nodes,core=qe_cores,name='qe_%s_proj'%prefix,walltime=walltime)
    qe.add_command('module load %s'%mod_espresso)
    qe.add_command('python proj_mote2.py -c -a')
    qe.run()


if args.analyse:
    pxml = ProjwfcXML(prefix,path=folder)
    print pxml
    print pxml.proj.dtype
    print pxml.proj.shape
    print "Writting projections"
    pxml.write_proj()
    print "done!"

if args.plot:
    pxml = ProjwfcXML(prefix,path=folder)
    print pxml
    no = 30
    l1 = list(xrange(no))
    l2 = list(xrange(no,2*no))
    l3 = list(xrange(2*no,3*no))

    s = [0,16]
    s = np.array(s)
    p = [1,2,3,17,18,19]
    p = np.array(p)

    sa = list(s)+list(s+no)+list(s+no*2)
    pa = list(p)+list(p+no)+list(p+no*2)
    #da = list(d)+list(d+no)+list(d+no*2)
    s = sa
    p = pa
    #d = da
    print s
    print p
    #print d

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(30,10))
    #for n,(orb,title) in enumerate(zip([l1,l2,l3],['layer1','layer2','layer3'])):
    for n,(orb,title) in enumerate(zip([s,p],['s','p'])):
        ax = plt.subplot(1,3,n+1)
        plt.title(title)

        #plot band lines
        eig = pxml.get_eigen()
        print eig
        ax.plot(eig-pxml.fermi,c='b',zorder=0)

        #plot character
        pxml.plot_eigen(ax,path=path,selected_orbitals=orb)

        ax.set_ylim([-7,4])
    plt.show()

