# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import numpy as np


def abs2(x):
    return x.real**2 + x.imag**2
 
def lorentzian(x,x0,g):
    height=1./(np.pi*g)
    return height*(g**2)/((x-x0)**2+g**2)

def gaussian(x,x0,s,max_exp=50.,min_exp=-100.):
    height=1./(np.sqrt(2.*np.pi)*s)
    argument=-0.5*((x-x0)/s)**2
    #Avoiding undeflow errors...
    np.place(argument,argument<min_exp,min_exp)
    return height*np.exp(argument)

def boltzman_f(Eb, Bose_Temp):
    kb = 8.61733326*10**-5
    return np.exp(-Eb/(kb*Bose_Temp))


