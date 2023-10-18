import sys
import numpy as np
from yambopy.units import ha2ev

def damp_it(ft, time, t_initial, damp_type="LORENTZIAN", damp_factor=0.1/ha2ev):
    ft_damped=np.empty_like(ft)

    if damp_type.upper() == "LORENTZIAN":
        ft_damped[:]=ft[:]*np.exp(-abs(time[:]-t_initial)*damp_factor)
    elif damp_type.upper() == "GAUSSIAN":
        ft_damped[:]=ft[:]*np.exp(-(time[:]-t_initial)**2*damp_factor**2)
    else:
        print("Wrong damping type ")
        sys.exit(0)

    return ft_damped
