#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC, AMS, FP, RR
#
# This file is part of the yambopy project
#
import numpy as np

def CmplxType(var):
    """ Distinguish between double and float for storing residuals and eigenvector with the same precision as in the Yambo database
    """
    if var.dtype=='float32':   return np.complex64
    elif var.dtype=='float64': return np.complex128
    else: raise TypeError('\n[ERROR] Variable type not recognized. It should be either float (float32) or double (float64).\n')
