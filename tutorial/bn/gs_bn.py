#
# Author: Henrique Pereira Coutada Miranda
# Run a Silicon groundstate calculation using Quantum Espresso
#
from __future__ import print_function, division
from pwpy.inputfile import *
from pwpy.outputxml import *
import argparse

#
# Create the input files
#
def get_inputfile():
    """ Define a Quantum espresso input file for boron nitride
    """ 
    qe = PwIn()
    qe.atoms = [['N',[0.0,0.0,0.5]],
                ['B',[1/3,2/3,0.5]]]
    qe.atypes = {'B': [10.811, "B.pbe-mt_fhi.UPF"],
                 'N': [14.0067,"N.pbe-mt_fhi.UPF"]}

    qe.control['prefix'] = "'bn'"
    qe.control['wf_collect'] = '.true.'
    qe.system['celldm(1)'] = 4.7
    qe.system['celldm(3)'] = 12/qe.system['celldm(1)']
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 2
    qe.system['ibrav'] = 4
    qe.kpoints = [12, 12, 1]
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
def scf():
    if not os.path.isdir('scf'):
        os.mkdir('scf')
    qe = get_inputfile()
    qe.control['calculation'] = "'scf'"
    qe.write('scf/bn.scf')
 
#nscf
def nscf(kpoints,folder):
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

def phonon(kpoints,qpoints,folder):
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
    e = EspressoXML('bn',path=pathin)
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
    parser.add_argument('-r' ,'--relax',       action="store_false", help='Don\'t run structural relaxation')
    parser.add_argument('-s' ,'--scf',         action="store_false", help='Don\'t run self-consistent calculation')
    parser.add_argument('-n' ,'--nscf',        action="store_false", help='Don\'t run non-self consistent calculation')
    parser.add_argument('-n2','--nscf_double', action="store_false", help='Don\'t run non-self consistent calculation for the double grid')
    parser.add_argument('-p' ,'--phonon',      action="store_false", help='Don\'t run phonon calculation')
    parser.add_argument('-t' ,'--nthreads',    action="store_true", help='Number of threads', default=2 )
    args = parser.parse_args()

    # create input files and folders
    relax()
    scf()
    nscf([6,6,6], 'nscf')
    nscf([12,12,1], 'nscf_double')
    phonon([12,12,1],[3,3,1], 'phonon')

    if args.relax:
        print("running relax:")
        os.system("cd relax; mpirun -np %d pw.x -inp bn.scf > relax.log"%args.nthreads)  #relax
        update_positions('relax','scf') 
        print("done!")

    if args.scf:
        print("running scf:")
        os.system("cd scf; mpirun -np %d pw.x -inp bn.scf > scf.log"%args.nthreads)  #scf
        print("done!")
   
    if args.nscf: 
        print("running nscf:")
        os.system("cp -r scf/bn.save nscf/") #nscf
        os.system("cd nscf; mpirun -np %d pw.x -inp bn.nscf > nscf.log"%args.nthreads) #nscf
        print("done!")

    if args.nscf_double: 
        print("running nscf_double:")
        os.system("cp -r scf/bn.save nscf_double/") #nscf
        os.system("cd nscf_double; mpirun -np %d pw.x -inp bn.nscf > nscf_double.log"%args.nthreads) #nscf
        print("done!")
    
    if args.phonon:
        print("running phonon:")
        os.system("cp -r scf/bn.save phonon/")
        os.system("cd phonon; mpirun -np %d ph.x -inp bn.ph > phonon.log"%args.nthreads) #phonon
        os.system("cd phonon; dynmat.x -inp bn.dynmat > dynmat.log") #matdyn
        print("done!")

