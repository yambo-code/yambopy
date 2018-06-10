# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
class YamboJson():
    """
    Read a json file produced by YamboOut
    """
    def __init__(self,filename):
        self.init = False

    @classmethod
    def from_file(filename):

        #open json file
        with open(filename,"r") as f:
            j = json.load(d)

        return YamboJson.from_dict(j)

    def from_dict(j):

        yj = YamboJson()
    
        yj.data =      j["data"]
        yj.tags =      j["tags"]
        yj.runtime =   j["runtime"]
        yj.inputfile = j["inputfile"]
        yj.lattice =   j["lattice"]
        yj.alat =      j["alat"]
        yj.kpts_iku =  j["kpts_iku"]
        yj.sym_car =   j["sym_car"]
        yj.atompos =   j["atompos"]
        yj.atomtype =  j["atomtype"]

        return yj
     
