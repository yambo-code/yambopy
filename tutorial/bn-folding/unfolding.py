#############################################################################################
# Unfolding of the electronic structure.
# Program adapted for reading values of Quantum Espresso.
#
# Authors: Alejandro Molina-Sanchez and Henrique Miranda
# Version of 24 of February of 2014.
#############################################################################################

from qepy import PwXML
from qepy import Unfolding
from numpy import array, dot, sqrt, cross
import matplotlib.pyplot as plt
from qepy import *

prefix_pc = 'pc'
prefix_sc = 'sc'
pc = PwXML(prefix=prefix_pc)
sc = PwXML(prefix=prefix_sc)

fold = Unfolding(prefix_pc=prefix_pc,path_pc='.',prefix_sc=prefix_sc,path_sc='.')

p = Path([ [[0.0, 0.0, 0.0],'G'],
           [[0.5, 0.0, 0.0],'M']], [11])

ax = plt.subplot(1,1,1)
fold.plot_eigen_ax(ax,path=p)
plt.show()
