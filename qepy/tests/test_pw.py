# Copyright (C) 2018 Henrique Pereira Coutada Miranda
# All rights reserved.
#
# This file is part of yambopy
#
import unittest
import os
from qepy.pw import PwIn
from yambopy.data.structures import Si

class TestPwIn(unittest.TestCase):

    def test_pwin(self):
        pwi = PwIn.from_structure_dict(Si) 
        pwi.set_nscf(10)
        print(pwi)

        pwi.cell_parameters = [[1,0,0],[0,1,0],[0,0,1]]
        pwi.ibrav = 0
        print(pwi)

if __name__ == '__main__':
    unittest.main()
