#
# Author: Henrique Pereira Coutada Miranda
# Run a GW calculation using Yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *
import argparse

#pack the files in .json files
pack_files_in_folder('gw_conv')

#plot the results using yambm analyser
ya = YamboAnalyser('gw_conv')
print(ya)
print('kpoints along a path')

path = [[0.5,   0,   0],
        [  0,   0,   0],
        [  0, 0.5, 0.5],
        [1.0, 1.0, 1.0]]
ya.get_path(path,'reference.json')

print('done!')
