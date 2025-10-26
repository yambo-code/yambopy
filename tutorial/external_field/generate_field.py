import sys
import argparse
from yambopy.units import ha2ev
"""
In this example we show how to generate and external field
in a format that readable from yambo_rt and yambo_nl.
We give two example SIN and QSSIN functions
"""

#
# An external field should provide f(t), f'(t), f''(t) and f(w)
# I put all paramters of the field inside their repecitve functions
# If you are not able to compute analitically f(w) you can use FFT
#

def sin_field(t_steps,w_steps):
    # In this example we generate a field with frequency 2.0 eV
    freq=2.0/ha2ev 
    f_t   =1.0/freq*cos(freq*t_steps)   # A(t) 
    fp_t  =sin(freq*t_steps)            # A'(t)
    fpp_t =-freq*sin(freq*t_steps)      # A''(t)
    f_w   =0
    return f_t,fp_t,fpp_t,f_w


def qsin_field(t_steps,w_steps):
    # In this example we generate a field with frequency 2.0 eV
    freq=2.0/ha2ev 
    f_t   =1.0/freq*cos(freq*t_steps)   # A(t) 
    fp_t  =sin(freq*t_steps)            # A'(t)
    fpp_t =-freq*sin(freq*t_steps)      # A''(t)
    f_w   =0
    return f_t,fp_t,fpp_t,f_w
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='External field generation',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--tstep',  type=float, default=0.01,      help='Time step (fs)')
    parser.add_argument('-t', '--trange', type=float, nargs=2, default=[0.0,80.0],help='Time range (fs)')
    parser.add_argument('-s', '--tstart', type=float, default=0.01,      help='Initial time (fs)')
    parser.add_argument('-v', '--versor', type=float, nargs=3, default=[1.0,0.0,0.0],help='Field versor')
    parser.add_argument('-f', '--fname',  type=str, default='SIN', help='Field name (SIN | QSIN) ')
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

