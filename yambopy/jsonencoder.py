#
# Copyright (C) 2017 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
#
"""import a special json encoder to do slightly readable files"""
import json
import re
import numpy as np

class YambopyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            if np.iscomplexobj(obj):
                return {"real":obj.real, "imag":obj.imag}
            else:
                return obj.tolist()
        return(json.JSONEncoder.default(self, obj))

def JsonDumpers(data):
    s = json.dumps(data,cls=YambopyEncoder,indent=1)
    return re.sub('\s+([0-9\]\-])','\\1',s)

def JsonDumper(data,filename):
    with open(filename,'w') as f:
        f.write(JsonDumpers(data))
