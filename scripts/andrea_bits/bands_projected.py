import numpy as np
from matplotlib import pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="file",required=True,help="out filename")
args = parser.parse_args()

data = np.loadtxt(args.file).T
#np.save('fort.100', data.T)
#data = np.load('fort.100.npy')
point_frac = 1
data_plot = data[:,::point_frac]

plt.figure(0,figsize=[5,10])
plt.scatter(data_plot[0], data_plot[1]-0.3, c=data_plot[2], edgecolors='none', cmap=plt.cm.get_cmap('PuBu'), vmax=max(data_plot[2])/2.)
plt.scatter(-data_plot[0], data_plot[1]-0.3, c=data_plot[2], edgecolors='none', cmap=plt.cm.get_cmap('PuBu'), vmax=max(data_plot[2])/2.)
#print dir(plt.gca())

plt.gca().set_ylabel('$E-E_F$ (eV)')
plt.yticks([0.0, -0.5, -1.0, -1.5, -2.0])
plt.ylim([-2.2,0.])
plt.gca().set_yticklabels(['0.0','-0.5', '-1.0', '-1.5', '-2.0'])

#plt.title('Armchair Direction')
#plt.gca().set_xlabel('$k_x$ $(\AA^{-1})$')
#pref = 0.529
#plt.xticks([-0.6*pref, -0.4*pref, -0.2*pref, 0.0, 0.2*pref, 0.4*pref,0.6*pref])
#plt.xlim([-0.71*pref,0.51*pref])
#plt.gca().set_xticklabels(['-0.6','-0.4','-0.2','0','0.2', '0.4', '0.6'])

plt.title('Zig-zag Direction')
plt.gca().set_xlabel('$k_x$ $(\AA^{-1})$')
pref = 0.529
plt.xticks([-0.6*pref, -0.4*pref, -0.2*pref, 0.0, 0.2*pref, 0.4*pref,0.6*pref])
plt.xlim([-0.71*pref,0.71*pref])
plt.gca().set_xticklabels(['-0.6','-0.4','-0.2','0','0.2', '0.4', '0.6'])

#plt.colorbar().remove
plt.savefig(args.file+'.pdf', format='pdf')
plt.savefig(args.file+'.png', format='png', dpi=300)
plt.show()
