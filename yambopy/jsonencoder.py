"""
import a special json encoder to do slightly readable files
"""
# Copyright (C) 2017 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
import json
import re
import numpy as np

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist() # or map(int, obj)
        return(json.JSONEncoder.default(self, obj))

def JsonDumpers(data):
    s = json.dumps(data,cls=MyEncoder,indent=1)
    return re.sub('\s+([0-9\]\-])','\\1',s)

def JsonDumperf(data,f):
    f.write(JsonDumpers(data))

def JsonDumper(data,filename):
    f = open(filename,'w')
    f.write(JsonDumpers(data))
    f.close()
