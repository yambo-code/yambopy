import sys
import argparse
import numpy as np
from yambopy.units import ha2ev,fs2aut
from yambopy.nl.compose_field   import Compose_Field
"""
In this example we show how to generate and external field
in a format that readable from yambo_rt and yambo_nl.
We give two example SIN and QSSIN functions
"""

#
# An external field should provide f(t), f'(t), f''(t) 
# I put all paramters of the field inside their repecitve functions
#

def sin_field(t_start,t_steps):
    # In this example we generate a field with frequency 2.0 eV
    freq=2.0/ha2ev 
    f_t   =-1.0/freq*(np.cos(freq*(t_start-t_steps))-1.0)   # A(t) 
    fp_t  = np.sin(freq*(t_start-t_steps))                  # A'(t)
    fpp_t = freq*np.cos(freq*(t_start-t_steps))             # A''(t)
    f_w   = 0
    return  f_t,fp_t,fpp_t 


#def qsin_field(t_start,t_steps,w_steps):
#    # In this example we generate a field with frequency 2.0 eV
#    freq=2.0/ha2ev ddd
#    f_t   =1.0/freq*cos(freq*(t_steps-t_start))   # A(t) 
#    fp_t  =sin(freq*(t_steps-t_start))            # A'(t)
#    fpp_t =-freq*sin(freq*(t_steps-t_start))      # A''(t)
#    f_w   =0
#    return f_t,fp_t,fpp_t,f_w
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='External field generation',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--tstep',  type=float, default=0.01,      help='Time step (fs)')
    parser.add_argument('-t', '--trange', type=float, nargs=2, default=[0.0,80.0],help='Time range (fs)')
    parser.add_argument('-s', '--tstart', type=float, default=0.01,      help='Initial time (fs)')
    parser.add_argument('-v', '--versor', type=float, nargs=3, default=[1.0,0.0,0.0],help='Field versor')
    parser.add_argument('-i', '--fint',   type=float, default=1000.0,      help='Field Intensity [kW/cm^2]')
    parser.add_argument('-f', '--fname',  type=str, help='Field name (SIN | QSIN) ')
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    if args.fname == None:
        print("You have to specify the field name ")
        sys.exit(1)

    t_start=args.tstart*fs2aut
    t_step =args.tstep*fs2aut
    t_range=np.arange(args.trange[0]*fs2aut,args.trange[1]*fs2aut,t_step)
    f_amp  =
    
    print("\n\n * * * Generate and external field for yambo_rt/yambo_nl * * * \n\n")
    print("Field name : ",args.fname)
    print("Time range : ",args.trange,"[fs]")
    print("Time step  : ",args.tstep, "[fs]")
    print("Start time : ",args.tstart,"[fs]")

    if args.fname == "SIN":
        a_pot=sin_field(t_start,t_range)
    elif args.fname == "QSIN":
        a_pot=qsin_field(t_start,t_range)

    field_fname="ext_field.txt"
    Compose_Field(a_pot,t_start,t_step,t_range,field_fname)

