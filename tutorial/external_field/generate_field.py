import sys
import argparse
"""
In this example we show how to generate and external field
in a format that readable from yambo_rt and yambo_nl.
We give two example SIN and QSSIN functions
"""

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

