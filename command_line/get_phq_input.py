import numpy as np
import argparse

"""
Script to update the explicit list of q-points in a ph input file (ldisp=.false., qplot=.true.).

- Reads the output of the scf calculation and the ph input.

Usage:
 :: -pw,--pwout='path/to/pw/output/file
 :: -ph,--phin='path/to/ph/input/file'

"""
def get_phq_input(pw_out,ph_in):

    # Read pw output
    with open(pw_out) as pw: lines = pw.readlines()
    Nlines = len(lines)

    for i in range(Nlines):
        if "number of k points=" in lines[i]:
            Nk = int( lines[i].split()[-1] )
            i_start = i+2
            break

    # Generate list in PH/yambo format
    kpts = []
    qpts = np.zeros([Nk,4])
    for ik in range(i_start,i_start+Nk): kpts.append(lines[ik].split("), wk")[0].split("= (")[1])
    for ik in range(Nk):
        for i in range(3): qpts[ik,i]=-float(kpts[ik].split()[i])

    # Check ph input and write
    with open(ph_in) as ph: lines = ph.readlines()

    if str(Nk) != lines[-1].strip():
        print("[ERROR] Check last lines of ph input")
        print("        - qpoints already there?")
        print("        - empty newline?")
        print("        - wrong Nk? (should be %d)"%Nk)
    else:
        ph = open(ph_in, 'a')
        for ik in range(Nk):
            qx, qy, qz = qpts[ik,0], qpts[ik,1], qpts[ik,2]
            ph.write('  %.7f  %.7f  %.7f  1\n'%(qx,qy,qz))
        ph.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Append explicit qpoints to ph input')
    parser.add_argument('-pw','--pwout', type=str, help='Path to pw (scf) output file',required=True)
    parser.add_argument('-ph','--phin', type=str,help='Path to ph (dvscf or elph) input file',required=True)
    args = parser.parse_args()

    pwout = args.pwout
    phin  = args.phin

    get_phq_input(pwout,phin)
