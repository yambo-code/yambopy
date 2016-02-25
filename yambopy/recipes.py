# Copyright (C) 2015 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
from yambopy.inputfile import YamboIn
import os

def breaking_symmetries(efield1,efield2=[0,0,0],folder='.',RTS=True):
# Breaks the symmetries for a given field.
# Second field used in circular polarized pump configuration
# RTS : Remove time symmetry is set True by default
  os.system('mkdir -p %s'%folder)
  os.system('cp -r database/SAVE %s'%folder)
  ypp = YamboIn('ypp_ph -n -V all',folder=folder,filename='ypp.in')
  ypp['Efield1'] = efield1 # Field in the X-direction
  ypp['Efield2'] = efield2 # Field in the X-direction
  if RTS:
    ypp.arguments.append('RmTimeRev')   # Remove Time Symmetry
  ypp.write('%s/ypp.in'%folder)
  os.system('cd %s ; ypp_ph -F ypp.in'%folder )
  os.system('cd %s ; cd FixSymm; yambo '%folder )
  os.system('rm -r %s/SAVE'%folder)
  os.system('mv %s/FixSymm/SAVE %s/'%(folder,folder))
  os.system('rm -r %s/FixSymm'%folder)
