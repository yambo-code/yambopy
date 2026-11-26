"""
Generate an external field in the format read by yambo_rt / yambo_nl
(field kind "FROM_FILE <file>" in the real-time input).

Usage
-----
    python3 generate_field.py -f GAUSS -d 0.005 -t 0 120 -s 0.0 -o EXTFIELD1_P1.time

Two example fields are implemented: SIN and GAUSS (see the functions
sin_field and gauss_field below). To create your own field write a new
function that returns the three arrays  a(t), a'(t), a''(t)  and add it
to the dictionary FIELDS.

Conventions
-----------
* a(t) is the vector potential, the electric field is  E(t) = -E_amp * a'(t).
* The field is saved WITHOUT amplitude: the amplitude is obtained by Yambo
  from the intensity (Field1_Int) of the real-time input.
* The polarization (versor) is NOT saved: it is set by Field1_Dir in the
  real-time input.
* The time inside the functions is in atomic units.

Format of the output file
-------------------------
    first row : N_steps  dt_f             (number of rows, time step in ATOMIC UNITS)
    next rows : t  a(t)  a'(t)  a''(t)    (t in FEMTOSECONDS, a, a', a'' in a.u.)

Time step
---------
The time step -d (dt_f) must be compatible with the one of the simulation (dt):
    EULER and INVINT                 dt_f = dt/n
    CRANKNIC, RK2EXP and RK2         dt_f = dt/(2n)
    RK4 and RK4EXP                   dt_f = dt/(6n)
A file with dt_f = dt/6 is valid for all the integrators.

Initial time
------------
The starting time of the field is given by Field1_Tstart in the real-time
input. If you also use -s here the shift is applied twice: use -s 0.0.
(The default of -s is 0.01 fs only for backward compatibility.)
"""

import argparse
import math
import numpy as np
from scipy import special
from yambopy.units import ha2ev, fs2aut

# --------------------------------------------------------------------------
# Parameters of the example fields (edit them or write your own function)
# --------------------------------------------------------------------------
SIN_FREQ_EV    = 2.0   # frequency of the SIN field [eV]
GAUSS_SIGMA_FS = 5.0   # sigma of the Gaussian [fs]; the pulse is centred at 3*sigma


def sin_field(t_steps):
    """SIN field with frequency SIN_FREQ_EV.
    t_steps in atomic units, returns a(t), a'(t), a''(t)."""
    freq  = SIN_FREQ_EV/ha2ev
    a_t   = -1.0/freq*(np.cos(freq*t_steps)-1.0)   # a(t)
    ap_t  = np.sin(freq*t_steps)                   # a'(t)
    app_t = freq*np.cos(freq*t_steps)              # a''(t)
    return a_t, ap_t, app_t


def gauss_field(t_steps):
    """Gaussian field: a'(t) is a Gaussian of width GAUSS_SIGMA_FS centred at 3*sigma.
    t_steps in atomic units, returns a(t), a'(t), a''(t)."""
    sigma = GAUSS_SIGMA_FS*fs2aut
    T_0   = 3.0*sigma
    Expf  = np.exp(-(t_steps-T_0)**2/(2.0*sigma**2))
    a_t   = sigma*np.sqrt(math.pi/2.0)*(special.erf((t_steps-T_0)/(sigma*np.sqrt(2.0)))+1.0)
    ap_t  = Expf
    app_t = -Expf*(t_steps-T_0)/sigma**2
    return a_t, ap_t, app_t


# name used with the option -f  ->  function that generates the field
FIELDS = {"SIN": sin_field, "GAUSS": gauss_field}


def write_field(file_name, t_steps, dt_au, a_pot):
    """Write the field on file, in the format described at the top of this script.
    t_steps and dt_au in atomic units (t is written in fs)."""
    data = np.column_stack((t_steps/fs2aut, a_pot[0], a_pot[1], a_pot[2]))
    with open(file_name, "w") as f:
        f.write(str(t_steps.size)+"     "+str(dt_au)+" \n")
        np.savetxt(f, data, fmt="%4.8e", delimiter="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='External field generation for yambo_rt/yambo_nl',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fname',  type=str, required=True, choices=sorted(FIELDS),
                        help='Field name')
    parser.add_argument('-d', '--tstep',  type=float, default=0.01,
                        help='Time step dt_f of the field (fs)')
    parser.add_argument('-t', '--trange', type=float, nargs=2, default=[0.0, 80.0],
                        metavar=('T_INIT', 'T_FINAL'), help='Time range (fs)')
    parser.add_argument('-s', '--tstart', type=float, default=0.01,
                        help='Shift of the field (fs). Use 0.0 if you set Field1_Tstart in the Yambo input')
    parser.add_argument('-o', '--fout',   type=str, default="EXTFIELD1_P1.time",
                        help='Name of the output file')
    parser.add_argument('-v', '--versor', type=float, nargs=3, default=[1.0, 0.0, 0.0],
                        help='IGNORED (kept for backward compatibility): the polarization is set by Field1_Dir in the Yambo input')
    args = parser.parse_args()

    t_start = args.tstart*fs2aut
    t_step  = args.tstep*fs2aut
    t_range = np.arange(args.trange[0]*fs2aut, args.trange[1]*fs2aut+t_step/2.0, t_step)

    print("\n * * * Generate an external field for yambo_rt/yambo_nl * * *\n")
    print("Field name : ", args.fname)
    print("Time range : ", args.trange, "[fs]")
    print("Time step  : ", args.tstep, "[fs]")
    print("Start time : ", args.tstart, "[fs]")
    print("Output file: ", args.fout)

    a_pot = FIELDS[args.fname](t_range-t_start)

    write_field(args.fout, t_range, t_step, a_pot)
