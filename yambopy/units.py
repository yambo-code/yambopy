import numpy as np

I = complex(0,1)
ha2ev  = 27.2113834 
ev2cm1 = 8065.5440044136285
bohr2ang = 0.52917720859
core_cnst =2.418884326505
fs2aut    =100.0/core_cnst 
as2aut    =0.1/core_cnst
ps2aut    =100000.0/core_cnst
speed_of_light=137.03599911
autime2s= 2.418884326505e-17
m_e = 0.510998950e+06 # electron mass in eV
KB = 8.617330337217213e-05/ha2ev # Ha/K
hbar = 6.582119569e-16 #eVs

SVCMm12VMm1 =29.98*np.power(10,3,dtype=np.double) #  from [statV/cm] to [V/m]
VMm12SVCMm1 =1.0/SVCMm12VMm1
AU2VMm1     =5.14220632*np.power(10,11,dtype=np.double)


atomic_mass = [   None,      1.00794,    4.002602,     6.941,   9.012182,
                10.811,      12.0107,     14.0067,   15.9994, 18.9984032,
               20.1797,  22.98976928,      24.305,26.9815386,    28.0855,
             30.973762,       32.065,      35.453,    39.948,    39.0983,
                40.078,    44.955912,      47.867,   50.9415,    51.9961,
             54.938045,       55.845,   58.933195,   58.6934,     63.546,
                 65.38,       69.723,       72.64,   74.9216,      78.96,
                79.904,       83.798,     85.4678,     87.62,   88.90585,
                91.224,     92.90638,       95.96,      None,     101.07,
              102.9055,       106.42,    107.8682,   112.411,    114.818,
                118.71,       121.76,       127.6, 126.90447,    131.293,
           132.9054519,      137.327,   138.90547,   140.116,  140.90765,
               144.242,         None,      150.36,   151.964,     157.25,
             158.92535,         162.5,  164.93032,   167.259,  168.93421,
               173.054,      174.9668,     178.49, 180.94788,     183.84,
               186.207,        190.23,    192.217,   195.084, 196.966569,
                200.59,      204.3833,      207.2,  208.9804,       None,
                  None,          None,       None,      None,       None,
             232.03806,     231.03588,  238.02891,      None,       None,
                  None,          None,       None,      None,       None,
                  None,          None,       None,      None,       None,
                  None,          None,       None,      None,       None,
                  None,          None,       None,      None,       None,
                  None,          None,       None,      None]

chemical_symbols = ['X',  'H',  'He', 'Li', 'Be',
                    'B',  'C',  'N',  'O',  'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si',
                    'P',  'S',  'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V',  'Cr',
                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se',
                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I',  'Xe',
                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
                    'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U',  'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']
