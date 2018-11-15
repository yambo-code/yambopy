# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from yamboparser import *
import os

class YamboHFDB():
    """
    Class to read yambo ndb.HF_and_locXC files
    """
    def __init__(self,filename='ndb.HF_and_locXC',folder='.'):
        """
        Read a QP file using the yamboparser
        """
        self.folder = folder
        self.filename = filename
        if os.path.isfile('%s/%s'%(folder,filename)):
            self.yfile = YamboFile(filename,folder)
        else:
            raise ValueError('File %s/%s not found'%(folder,filename))

        qps = self.yfile.data
        nkpoints = qps['nkpoints']
        nbands = qps['nbands']
        self.sx = qps['Sx'].reshape(nkpoints,nbands)
        self.vxc = qps['Vxc'].reshape(nkpoints,nbands)

    def __str__(self):
        s = ""
        s += "nqps:     %d\n"%self.nqps
        s += "nkpoints: %d\n"%self.nkpoints
        s += "nbands:   %d\n"%self.nbands
        return s

