# Copyright (C) 2018 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
#
#
import os

def get_pseudo_path(pseudo_name,absolute_paths):
    """Check if a pseudo with pseudo_name exists in any of the paths"""
    for path in absolute_paths:
        for filename in os.listdir(path):
            if pseudo_name == filename:
                return os.path.join(path,filename)
    return None
