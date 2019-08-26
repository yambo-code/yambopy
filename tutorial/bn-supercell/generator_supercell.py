from qepy import *
from yambopy.lattice import *
import argparse
import os
#test supercell class
filename1 = 'unit-cell.scf'
filename2 = 'unit-cell.scf'
modes_file, N_M = 'phonon.modes', 12 #*.modes file from PH phonon calculation at specific Q-point
#
# (A) python gen_supercell.py -d
# (B) python gen_supercell.py -nd
# (C) python gen_supercell.py -nd -disp
parser = argparse.ArgumentParser(description='Script to run supercell test')
parser.add_argument('-d' ,'--diagonal', action="store_true", help='Generate diagonal supercell')
parser.add_argument('-nd' ,'--nondiagonal', action="store_true", help='Generate non-diagonal supercell')
parser.add_argument('-disp','--displace', action="store_true", help='Displace supercell')
parser.add_argument('-r', '--remove', action="store_true", help='Remove previously generated files')
args = parser.parse_args()
# To clean the folder: python gen_supercell.py -r 
if args.remove: os.system('rm -f n_* d_* *_expanded')
#
#Diagonal
#
if args.diagonal:
    qe1 = PwIn.from_file(filename1)               #PwIn class
    R=[3,3,1]                           #[INPUT] Repetitions of unit cell in the lattice directions
    ysup_diag = Supercell(qe1)          #Supercell class
    ysup_diag.d_sup(R)                  #Generate diagonal supercell
    # Print supercell q-e input file
    # qe_d is the new supercell in PwIn() format
    calculation = ''.join(ysup_diag.qe_d.control['calculation'].split('\''))
    prefix = ''.join(ysup_diag.qe_d.control['prefix'].split('\''))
    ysup_diag.qe_d.write('d_%s.%s'%(prefix,calculation))
    print('supercell input file written.')
#
#Nondiagonal
#
if args.nondiagonal:
    qe2 = PwIn.from_file(filename2)           #PwIn class
    qe2.kpoints = [12,12,4]         #[OPTIONAL] I manually change the original k-point mesh of the input qe file (must be consistent with denominators Q[1])
    Q=[[1,-1,0],[3,6,2]]            #[INPUT] q-point to be folded at Gamma in fractional coord. in the BZ (Q[0]: numerators, Q[1]: denominators)
    ysup_ndiag = Supercell(qe2)     #Supercell class
    nd_sup = ysup_ndiag.nd_sup(Q)   #Generate non-diagonal supercell
    # Print supercell q-e input file
    # qe_nd is the new supercell in PwIn() format
    calculation = ''.join(ysup_ndiag.qe_nd.control['calculation'].split('\''))
    prefix = ''.join(ysup_ndiag.qe_nd.control['prefix'].split('\''))
    ysup_ndiag.qe_nd.write('n_%s.%s'%(prefix,calculation))
    print('supercell input file written.')
    # Displace nondiagonal
    if args.displace:
        ysup_ndiag.displace(modes_file,nd_sup,Temp=0.1) #Displace atoms: intensity is Temp [BOHR], sign and direction (standing wave with Q) given by modes_file
        # Print supercell q-e input files
        # modes_qe is the list of displaced supercell in PwIn() format along each phonon mode
        for mode in range(N_M): ysup_ndiag.modes_qe[mode].write('n_%s_mode%d.%s'%(prefix,mode+1,calculation))
        print('displaced supercell input files written.')

