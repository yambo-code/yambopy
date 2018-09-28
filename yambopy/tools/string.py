# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
def marquee(string,width=60,mark='='):
    return string.center(width).replace(" ",mark)
