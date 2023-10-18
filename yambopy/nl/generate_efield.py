import sys
import numpy as np
import scipy as special
import math  

def get_Efield_w(freqs,efield):
    if efield["name"] == "DELTA":
        efield_w=efield["amplitude"]*np.exp(1j*freqs[:]*efield["initial_time"])
    else:
        print("Fields different from Delta function not implemented yet")
        sys.exit(0)

    return efield_w

    
def Divide_by_the_Field(efield,order):
    
    
    if efield['name']=='SIN' or efield['name']=='SOFTSIN':
        if order !=0:
            divide_by_field=-2.*1j/efield['amplitude']**order
        elif order==0:
            divide_by_field=4/efield['amplitude']**2

    elif efield['name'] == 'QSSIN':
        # Approximate relations/does not work yet
        sigma=efied['width']
        T_0=10*sigma
        W_0=efield['frequency'][0]
        T = 2*np.pi/W_0
        E_w= math.sqrt(np.pi/2)*sigma*np.exp(-1j*W_0*T_0)*(special.erf((T-T_0)/math.sqrt(2.0)/sigma)+special.erf(T_0/sqrt(2.0)/sigma))
        
        if order!=0:
            divide_by_field = (2*1j/(E_w*efield['amplitude']))**order
        elif order==0:
            divide_by_field = 4*1j/(E_w*efield['amplitude']*np.conj(E_w))
    else:
        print('Electric field not implemented in Divide_by_the_Field!')
        sys.exit(0)

    return divide_by_field


