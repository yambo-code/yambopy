# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from yambopy import *
import os

def pack_files_in_folder(folder,save_folder=None,mask='',verbose=True):
    """
     Pack the output files in a folder to json files
    """
    if not save_folder: save_folder = folder
    print mask
    #pack the files in .json files
    for dirpath,dirnames,filenames in os.walk(folder):
        #check if the folder fits the mask
        print 'if %s in %s' %(mask, dirpath)
        if mask in dirpath:
            print 'OK'
            #check if there are some output files in the folder
            if ([ f for f in filenames if 'o-' in f ]):
                print dirpath
                y = YamboOut(dirpath,save_folder=save_folder)
                y.pack()

def breaking_symmetries(efield1,efield2=[0,0,0],folder='.',RmTimeRev=True):
# Breaks the symmetries for a given field.
# Second field used in circular polarized pump configuration
# RmTimeRev : Remove time symmetry is set True by default
  os.system('mkdir -p %s'%folder)
  os.system('cp -r database/SAVE %s'%folder)
  os.system('cd %s; yambo'%folder)
  ypp = YamboIn('ypp_ph -y -V all',folder=folder,filename='ypp.in')
  ypp['Efield1'] = efield1 # Field in the X-direction
  ypp['Efield2'] = efield2 # Field in the X-direction
  if RmTimeRev:
    ypp.arguments.append('RmTimeRev')   # Remove Time Symmetry
  ypp.write('%s/ypp.in'%folder)
  os.system('cd %s ; ypp_ph -F ypp.in'%folder )
  os.system('cd %s ; cd FixSymm; yambo '%folder )
  os.system('rm -r %s/SAVE'%folder)
  os.system('mv %s/FixSymm/SAVE %s/'%(folder,folder))
  os.system('rm -r %s/FixSymm'%folder)
