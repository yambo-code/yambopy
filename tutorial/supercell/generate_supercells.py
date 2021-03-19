from qepy import *
#from yambopy.lattice import *
import sys
import argparse
"""
In this example we construct supercells (sc) of bulk hBN. 
Starting from a unit cell (uc) input that is read from file, we build:

- diagonal supercell
- non-diagonal supercell folding a specific q-point Q
- displaced non-diagonal supercells along phonon eigenmodes at Q (read from file)
"""

def generate_diagonal_supercell(uc,R):
    """
    First case: diagonal supercell of size R
    """
    qe = PwIn.from_file(uc) #read the uc input

    sc = Supercell(qe) # initialize class
    sc.d_sup(R) # Generate supercell as PwIn object called 'qe_d'

    #name of output file
    suffix = ''.join(sc.qe_d.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_d.control['prefix'].split('\''))

    #write supercell to file
    sc.qe_d.write('sc_diagonal_%s.%s'%(prefix,suffix))

    print('Diagonal supercell written to file.')

def generate_nondiagonal_supercell(uc,Q,kpoints=None):
    """
    Second case: nondiagonal supercell folding point Q
    """
    qe = PwIn.from_file(uc) #read the uc input
    if kpoints is not None: qe.kpoints = kpoints #Optional: manually change the original kpt mesh to ensure consistency    

    sc = Supercell(qe) # initialize class
    sc.nd_sup(Q) # Generate supercell as PwIn object called 'qe_nd'

    #name of output file
    suffix = ''.join(sc.qe_nd.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_nd.control['prefix'].split('\''))

    #write supercell to file
    sc.qe_nd.write('sc_nondiagonal_%s.%s'%(prefix,suffix))

    print('Nondiagonal supercell written to file.')   

def generate_displaced_supercells(uc,Q,eivs,kpoints=None):
    """
    Third case: displaced supercells along phonon modes at Q
    """
    qe = PwIn.from_file(uc) #read the uc input
    if kpoints is not None: qe.kpoints = kpoints #Optional: manually change the original kpt mesh to ensure consistency

    sc = Supercell(qe) # initialize class
    nd_atom_positions = sc.nd_sup(Q) # Generate supercell as PwIn object called 'qe_nd' getting new atomic positions

    # Displace atoms.
    # Intensity is Temp (in bohr)
    # Sign and direction (standing wave at Q) given by modes_file
    sc.displace(eivs,nd_atom_positions,Temp=0.1) # Generate list of displaced supercells as PwIn objects called 'modes_qe'
    N_modes = len(sc.modes_qe)

    #name of output file
    suffix = ''.join(sc.qe_nd.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_nd.control['prefix'].split('\''))

    #write supercells to file for each mode
    for mode in range(N_modes): sc.modes_qe[mode].write('sc_displaced_%s_mode%d.%s'%(prefix,mode+1,suffix)) 

    print('Displaced supercells written to file.')

def generate_displaced_unitcell(uc,eivs):
    """
    Fourth case: displaced cell at Q=0
    """
    qe = PwIn.from_file(uc) #read the uc input

    sc = Supercell(qe) # initialize class
    atom_positions = sc.d_sup([1,1,1]) # Generate "supercell" with size 1 as PwIn object called 'qe_d' getting atomic positions
    
    # Displace atoms.
    # Intensity is Temp (in bohr)
    # Sign and direction (standing wave at Q) given by modes_file
    sc.displace(eivs,atom_positions,Temp=0.1) # Generate list of displaced supercells as PwIn objects called 'modes_qe'
    N_modes = len(sc.modes_qe)

    #name of output file
    suffix = ''.join(sc.qe_d.control['calculation'].split('\''))
    prefix = ''.join(sc.qe_d.control['prefix'].split('\''))

    #write supercells to file for each mode
    for mode in range(N_modes): sc.modes_qe[mode].write('uc_GAMMA_displaced_%s_mode%d.%s'%(prefix,mode+4,suffix))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Supercell generation')
    parser.add_argument('-D', '--diagonal', action="store_true", help='Build diagonal supercell')
    parser.add_argument('-N', '--nondiagonal', action="store_true", help='Build nondiagonal supercell')
    parser.add_argument('-S', '--displaced_nondiagonal', action="store_true", help='Build displaced nondiagonal supercells')
    parser.add_argument('-G', '--displaced_gamma', action="store_true", help='Build displaced unit cell')
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    uc_filnm = 'uc.scf' #pw input
    R = [3,3,1] #diagonal supercell dimensions
    Q = [[1,-1,0],[3,6,2]] #qpoint to fold in nondiagonal supercell: Q=(1/3,-1/6,0)
    kpoints = [12,12,4]    #NB: these fractional crystal coordinates must exactly divide the kpoint mesh!
    modes_file = 'Q.modes' #file with phonon eigenmodes at Q. NB2: we use those already scaled with the atomic masses

    if args.diagonal: generate_diagonal_supercell(uc_filnm,R)

    if args.nondiagonal: generate_nondiagonal_supercell(uc_filnm,Q,kpoints)

    #NB3: If computing the eigenmodes at Q with qe, remember to use Q in *Cartesian* coordinates
    #     in both qe input and output!!
    #     For example, Q_crystal=(1/3,-1/6,0) ====> Q_cartesian = (1/3,0,0) in the BN hexagonal lattice

    if args.displaced_nondiagonal: generate_displaced_supercells(uc_filnm,Q,modes_file,kpoints)
    
    #NB4: For the following we are using the same modes_file as in the case above for simplicity, but it is *wrong*:
    #     Phonon modes at Gamma must be used in this case.
    if args.displaced_gamma: generate_displaced_unitcell(uc_filnm,modes_file)
