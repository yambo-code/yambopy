import sys
import argparse
import numpy as np
import math
from scipy import special
from yambopy.units import ha2ev,fs2aut,AU2KWCMm2

"""
In this example we show how to generate and external field
in a format that readable from yambo_rt and yambo_nl.
We give two example SIN and GAUSS functions

An external field should provide f(t), f'(t), f''(t) 
I put all paramters of the field inside their repecitve functions

*** IMPORTANT ***
The field is save without field intensity(amplitude) this is read
from the yambo_nl/yambo_rt input
"""

def sin_field(t_steps):
    # In this example we generate a field with frequency 2.0 eV
    freq=2.0/ha2ev 
    f_t   = -1.0/freq*(np.cos(freq*t_steps)-1.0)   # A(t) 
    fp_t  = np.sin(freq*t_steps)                  # A'(t)
    fpp_t = freq*np.cos(freq*t_steps)             # A''(t)
    return  f_t,fp_t,fpp_t 


def gauss_field(t_steps):
    # In this example we generate a guassian external field
    sigma=5.0*fs2aut # Field width
    T_0  =3.0*sigma  
    Expf =np.exp(-(t_steps-T_0)**2/(2.0*sigma**2) )
    f_t  =sigma*np.sqrt(math.pi/2.0)*(special.erf((t_steps-T_0)/(sigma*np.sqrt(2.0)))+1.0)
    fp_t =Expf
    fpp_t=-Expf*(t_steps-T_0)/sigma**2
    return  f_t,fp_t,fpp_t 

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='External field generation',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--tstep',  type=float, default=0.01,      help='Time step (fs)')
    parser.add_argument('-t', '--trange', type=float, nargs=2, default=[0.0,80.0],help='Time range (fs)')
    parser.add_argument('-s', '--tstart', type=float, default=0.01,      help='Initial time (fs)')
    parser.add_argument('-v', '--versor', type=float, nargs=3, default=[1.0,0.0,0.0],help='Field versor')
    parser.add_argument('-f', '--fname',  type=str, help='Field name (SIN | GAUSS) ')
    parser.add_argument('-o', '--fout',  type=str, default="EXTFIELD1_P1.time",help='External field file name')
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    if args.fname == None:
        print("You have to specify the field name ")
        sys.exit(1)

    t_start=args.tstart*fs2aut
    t_step =args.tstep*fs2aut
    t_range=np.arange(args.trange[0]*fs2aut,args.trange[1]*fs2aut+t_step/2.0,t_step)
    
    print("\n\n * * * Generate and external field for yambo_rt/yambo_nl * * * \n\n")
    print("Field name : ",args.fname)
    print("Time range : ",args.trange,"[fs]")
    print("Time step  : ",args.tstep, "[fs]")
    print("Start time : ",args.tstart,"[fs]")

    if args.fname == "SIN":
        a_pot=sin_field(t_range-t_start)
    elif args.fname == "GAUSS":
        a_pot=gauss_field(t_range-t_start)
    else:
        print("Unknown external field ")
        sys.exit(1)

    data = np.column_stack((t_range/fs2aut,a_pot[0],a_pot[1],a_pot[2]))
    with open(args.fout, "w") as f:
        f.write(str(t_range.size)+"     "+str(t_step)+" \n")
        np.savetxt(f, data, fmt="%4.8e", delimiter="\t")
