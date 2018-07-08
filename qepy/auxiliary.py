# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina-Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
from builtins import range
from builtins import object
import numpy as np
import os

def float_from_string(x):
  """
  Convert a string in a float 
  """
  y=[]
  for t in x.split():
    try:
      y.append(float(t))
    except ValueError:
      pass
  return y


