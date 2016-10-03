#
# Author: Henrique Pereira Coutada Miranda
# Run a Silicon groundstate calculation using Quantum Espresso
#
from __future__ import print_function, division
import sys
from qepy import *
import argparse

kpoints = [24,24,1]
kpoints_double = [18,18,1]
qpoints = [3,3,1]
pw = 'pw.x'

# create the input files
def get_inputfile():
    """ Define a Quantum espresso input file for boron nitride
    """ 
    qe = PwIn()
    qe.atoms = [['N',[0.0,0.0,0.5]],
                ['B',[1/3,2/3,0.5]]]
    qe.atypes = {'B': [10.811, "B.pbe-mt_fhi.UPF"],
                 'N': [14.0067,"N.pbe-mt_fhi.UPF"]}

    qe.control['prefix'] = "'bn'"
    qe.control['verbosity'] = "'high'"
    qe.control['wf_collect'] = '.true.'
    qe.system['celldm(1)'] = 4.7
    qe.system['celldm(3)'] = 12/qe.system['celldm(1)']
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 2
    qe.system['ibrav'] = 4
    qe.kpoints = [9, 9, 1]
    qe.electrons['conv_thr'] = 1e-8
    return qe

#relax
def relax():
    if not os.path.isdir('relax'):
        os.mkdir('relax')
    qe = get_inputfile()
    qe.control['calculation'] = "'vc-relax'"
    qe.ions['ion_dynamics']  = "'bfgs'"
    qe.cell['cell_dynamics']  = "'bfgs'"
    qe.cell['cell_dofree']  = "'2Dxy'"
    qe.write('relax/bn.scf')

#scf
def scf(folder='scf'):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    qe = get_inputfile()
    qe.control['calculation'] = "'scf'"
    qe.write('%s/bn.scf'%folder)
 
#nscf
def nscf(kpoints,folder='nscf'):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    qe = get_inputfile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 30
    qe.system['force_symmorphic'] = ".true."
    qe.kpoints = kpoints
    qe.write('%s/bn.nscf'%folder)

def phonon(kpoints,qpoints,folder='phonon'):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    ph = PhIn()
    ph['nq1'],ph['nq2'],ph['nq3'] = qpoints
    ph['tr2_ph'] = 1e-12
    ph['prefix'] = "'bn'"
    ph['epsil'] = ".false."
    ph['trans'] = ".true."
    ph['fildyn'] = "'bn.dyn'"
    ph['fildrho'] = "'bn.drho'"
    ph['ldisp'] = ".true."
    ph.write('%s/bn.ph'%folder)

    md = DynmatIn()
    md['asr'] = "'simple'"
    md['fildyn'] = "'bn.dyn1'"
    md['filout'] = "'bn.modes'"
    md.write('%s/bn.dynmat'%folder)

def update_positions(pathin,pathout):
    """ update the positions of the atoms in the scf file using the output of the relaxation loop
    """
    e = PwXML('bn',path=pathin)
    pos = e.get_scaled_positions()

    q = PwIn('%s/bn.scf'%pathin)
    print("old celldm(1)", q.system['celldm(1)'])
    q.system['celldm(1)'] = e.cell[0][0]
    print("new celldm(1)", q.system['celldm(1)'])
    q.atoms = zip([a[0] for a in q.atoms],pos)
    q.write('%s/bn.scf'%pathout)

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r' ,'--relax',       action="store_true", help='Structural relaxation')
    parser.add_argument('-s' ,'--scf',         action="store_true", help='Self-consistent calculation')
    parser.add_argument('-n' ,'--nscf',        action="store_true", help='Non-self consistent calculation')
    parser.add_argument('-n2','--nscf_double', action="store_true", help='Non-self consistent calculation for the double grid')
    parser.add_argument('-p' ,'--phonon',      action="store_true", help='Phonon calculation')
    parser.add_argument('-t' ,'--nthreads',                         help='Number of threads', default=2 )
    args = parser.parse_args()
    nthreads = int(args.nthreads)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # create input files and folders
    relax()
    scf()
    nscf(kpoints)
    nscf(kpoints_double, folder='nscf_double')
    phonon(kpoints,qpoints)

    if args.relax:
        print("running relax:")
        os.system("cd relax; mpirun -np %d %s -inp bn.scf > relax.log"%(nthreads,pw))  #relax
        update_positions('relax','scf') 
        print("done!")

    if args.scf:
        print("running scf:")
        os.system("cd scf; mpirun -np %d %s -inp bn.scf > scf.log"%(nthreads,pw))  #scf
        print("done!")
   
    if args.nscf: 
        print("running nscf:")
        os.system("cp -r scf/bn.save nscf/") #nscf
        os.system("cd nscf; mpirun -np %d %s -nk %d -inp bn.nscf > nscf.log"%(nthreads,pw,nthreads)) #nscf
        print("done!")

    if args.nscf_double: 
        print("running nscf_double:")
        os.system("cp -r scf/bn.save nscf_double/") #nscf
        os.system("cd nscf_double; mpirun -np %d %s -inp bn.nscf > nscf_double.log"%(nthreads,pw)) #nscf
        print("done!")
    
    if args.phonon:
        print("running phonon:")
        os.system("cp -r scf/bn.save phonon/")
        os.system("cd phonon; mpirun -np %d %s -inp bn.ph > phonon.log"%(nthreads,ph)) #phonon
        os.system("cd phonon; dynmat.x -inp bn.dynmat > dynmat.log") #matdyn
        print("done!")

