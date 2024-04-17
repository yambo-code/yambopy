#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: RR
#
# This file is part of the yambopy project
#
import re
import numpy as np
from yambopy.lattice import replicate_red_kmesh, calculate_distances, get_path, car_red
from yambopy.dbs.latticedb import *
from yambopy.io.xsffile import *

def interlayer_power(xsf_hole_bottom_file, xsf_hole_top_file, zthreshold, c , fractional = False):
    xsf_hole_bottom = YamboXsf.read_xsf(xsf_hole_bottom_file)
    xsf_hole_top = YamboXsf.read_xsf(xsf_hole_top_file)        
    b_hole_b_contr, t_hole_b_contr = xsf_hole_bottom.contribution_twolayers(zthreshold,c, fractional )
    b_hole_t_contr, t_hole_t_contr = xsf_hole_top.contribution_twolayers(zthreshold,c, fractional )
    print('ciao',t_hole_b_contr+b_hole_b_contr+b_hole_t_contr+t_hole_t_contr)
    return ((t_hole_b_contr+b_hole_t_contr)/(b_hole_b_contr+t_hole_b_contr+t_hole_t_contr+b_hole_t_contr))
