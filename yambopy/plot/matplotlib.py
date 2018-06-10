# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
"""Check if we are on a headless system and in that case use Agg interface"""
import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    mpl.use('Agg')
from matplotlib import pyplot as plt
