# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#

"""
This file contains examples of structures
"""

__all__ = [
    "BN",
    "MoS2",
    "Si"
]

#BN
a = 4.7
c = 12
lattice = dict(ibrav=4,celldm1=4.7,celldm3=c/a)
atypes = dict(B=[10.811, "B.pbe-mt_fhi.UPF"],
              N=[14.0067,"N.pbe-mt_fhi.UPF"])

atoms = [['N',[ 0.0, 0.0,0.5]],
         ['B',[1./3,2./3,0.5]]]
BN = dict(lattice=lattice,atypes=atypes,atoms=atoms)

#MoS2
a = 5.838
c = 18
lattice = dict(ibrav=4,celldm1=5.838,celldm3=c/a) 
atypes = dict(Mo=[95.94, "Mo.pz-mt_fhi.UPF"],
              S =[32.065, "S.pz-mt_fhi.UPF"]) 
atoms = [['Mo',[2./3,1./3,          0.0]],
         [ 'S',[1./3,2./3, 2.92781466/c]],
         [ 'S',[1./3,2./3,-2.92781466/c]]]
MoS2 = dict(lattice=lattice,atypes=atypes,atoms=atoms)

#Si
lattice = dict(ibrav=2,celldm1=10.3)
atypes = dict(Si=[28.086,"Si.pbe-mt_fhi.UPF"])
atoms = [['Si',[0.125,0.125,0.125]],
         ['Si',[-.125,-.125,-.125]]]
Si = dict(lattice=lattice,atypes=atypes,atoms=atoms)

