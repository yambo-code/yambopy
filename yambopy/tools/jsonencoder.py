#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: HPC
#
# This file is part of the yambopy project
#
"""import a special json encoder to do slightly readable files"""
import json
import re
import numpy as np

class YambopyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.ndarray,np.number)):
            if np.iscomplexobj(obj):
                return [obj.real, obj.imag]
            else:
                return obj.tolist()
        return(json.JSONEncoder.default(self, obj))

def JsonLoaders(string):
    """ load json from string """
    return json.loads(filename)

def JsonLoader(filename):
    """ load json from file """
    with open(filename,'r') as f:
        data = json.load(f)
    return data

def JsonDumpers(data):
    """ dump dicitonary as string """
    s = json.dumps(data,cls=YambopyEncoder,indent=2)
    return re.sub('\s+([0-9\]\-])','\\1',s)

def JsonDumper(data,filename):
    """ dump dictionary as file """
    with open(filename,'w') as f:
        f.write(JsonDumpers(data))
