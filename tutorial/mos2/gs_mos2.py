#
# Author: Henrique Pereira Coutada Miranda
# Run a MoS2 groundstate calculation using Quantum Espresso
#
from __future__ import print_function, division
from pwpy.inputfile import *
from pwpy.outputxml import *
import argparse

#
# Create the input files
#
def get_inputfile():
    """ Define a Quantum espresso input file for MoS2
    """ 
    qe = PwIn()
    a = 5.83803209416
    c = 20
    qe.atoms = [['Mo',[2/3,1/3,0.0]],
                [ 'S',[1/3,2/3, 2.92781466/c]],
                [ 'S',[1/3,2/3,-2.92781466/c]]]
    qe.atypes = {'Mo': [10.811, "Mo.pz-mt_fhi.UPF"],
                 'S':  [14.0067, "S.pz-mt_fhi.UPF"]}

    qe.control['prefix'] = "'mos2'"
    qe.control['wf_collect'] = '.true.'
    qe.system['celldm(1)'] = a
    qe.system['celldm(3)'] = c/qe.system['celldm(1)']
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 3
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
    qe.write('relax/mos2.scf')

#scf
def scf():
    if not os.path.isdir('scf'):
        os.mkdir('scf')
    qe = get_inputfile()
    qe.control['calculation'] = "'scf'"
    qe.write('scf/mos2.scf')
 
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
    qe.write('%s/mos2.nscf'%folder)

def update_positions(pathin,pathout):
    """ update the positions of the atoms in the scf file using the output of the relaxation loop
    """
    e = EspressoXML('mos2',path=pathin)
    pos = e.get_scaled_positions()

    q = PwIn('%s/mos2.scf'%pathin)
    print("old celldm(1)", q.system['celldm(1)'])
    q.system['celldm(1)'] = e.cell[0][0]
    print("new celldm(1)", q.system['celldm(1)'])
    q.atoms = zip([a[0] for a in q.atoms],pos)
    q.write('%s/mos2.scf'%pathout)

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r' ,'--relax',       action="store_false", help='Don\'t run structural relaxation')
    parser.add_argument('-s' ,'--scf',         action="store_false", help='Don\'t run self-consistent calculation')
    parser.add_argument('-n' ,'--nscf',        action="store_false", help='Don\'t run non-self consistent calculation')
    parser.add_argument('-n2','--nscf_double', action="store_false", help='Don\'t run non-self consistent calculation for the double grid')
    parser.add_argument('-t' ,'--nthreads',    action="store_true", help='Number of threads', default=2 )
    args = parser.parse_args()

    # create input files and folders
    relax()
    scf()
    nscf([12,12,1], 'nscf')
    nscf([24,24,1], 'nscf_double')

    if args.relax:
        print("running relax:")
        os.system("cd relax; mpirun -np %d pw.x -inp mos2.scf > relax.log"%args.nthreads)  #relax
        update_positions('relax','scf') 
        print("done!")

    if args.scf:
        print("running scf:")
        os.system("cd scf; mpirun -np %d pw.x -inp mos2.scf > scf.log"%args.nthreads)  #scf
        print("done!")
   
    if args.nscf: 
        print("running nscf:")
        os.system("cp -r scf/mos2.save nscf/") #nscf
        os.system("cd nscf; mpirun -np %d pw.x -inp mos2.nscf > nscf.log"%args.nthreads) #nscf
        print("done!")

    if args.nscf_double: 
        print("running nscf_double:")
        os.system("cp -r scf/mos2.save nscf_double/") #nscf
        os.system("cd nscf_double; mpirun -np %d pw.x -inp mos2.nscf > nscf_double.log"%args.nthreads) #nscf
        print("done!")

