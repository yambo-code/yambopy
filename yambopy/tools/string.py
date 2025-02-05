#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC
#
# This file is part of the yambopy project
#
def marquee(string,width=60,mark='='):
    return string.center(width).replace(" ",mark)
