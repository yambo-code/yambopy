#
# Author: Henrique Pereira Coutada Miranda
# Run a Silicon groundstate calculation using Quantum Espresso
#
from __future__ import print_function
import sys
from qepy import *
import argparse

scf_kpoints  = [2,2,2]
nscf_kpoints = [3,3,3]
prefix = 'si'

p = Path([ [[1.0,1.0,1.0],'G'],
           [[0.0,0.5,0.5],'X'],
           [[0.0,0.0,0.0],'G'],
           [[0.5,0.0,0.0],'L']], [20,20,20])

#
# Create the input files
#
def get_inputfile():
    """ Define a Quantum espresso input file for silicon
    """
    qe = PwIn()
    qe.atoms = [['Si',[0.125,0.125,0.125]],
                ['Si',[-.125,-.125,-.125]]]
    qe.atypes = {'Si': [28.086,"Si.pbe-mt_fhi.UPF"]}

    qe.control['prefix'] = "'%s'"%prefix
    qe.control['wf_collect'] = '.true.'
    qe.system['celldm(1)'] = 10.3
    qe.system['ecutwfc'] = 30
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 1
    qe.system['ibrav'] = 2
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
    qe.kpoints = scf_kpoints
    qe.write('relax/%s.scf'%prefix)

#scf
def scf():
    if not os.path.isdir('scf'):
        os.mkdir('scf')
    qe = get_inputfile()
    qe.control['calculation'] = "'scf'"
    qe.kpoints = scf_kpoints
    qe.write('scf/%s.scf'%prefix)

#nscf
def nscf():
    if not os.path.isdir('nscf'):
        os.mkdir('nscf')
    qe = get_inputfile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 20
    qe.system['force_symmorphic'] = ".true."
    qe.kpoints = nscf_kpoints
    qe.write('nscf/%s.nscf'%prefix)

def bands():
    if not os.path.isdir('bands'):
        os.mkdir('bands')
    qe = get_inputfile()
    qe.control['calculation'] = "'bands'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 8
    qe.system['force_symmorphic'] = ".true."
    qe.ktype = 'crystal'
    qe.set_path(p)
    qe.write('bands/%s.bands'%prefix)

def update_positions(pathin,pathout):
    """ update the positions of the atoms in the scf file using the output of the relaxation loop
    """
    e = PwXML(prefix,path=pathin)
    pos = e.get_scaled_positions()

    q = PwIn('%s/%s.scf'%(pathin,prefix))
    print("old celldm(1)", q.system['celldm(1)'])
    q.system['celldm(1)'] = e.cell[0][2]*2
    print("new celldm(1)", q.system['celldm(1)'])
    q.atoms = zip([a[0] for a in q.atoms],pos)
    q.write('%s/%s.scf'%(pathout,prefix))

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r' ,'--relax',       action="store_true", help='Structural relaxation')
    parser.add_argument('-s' ,'--scf',         action="store_true", help='Self-consistent calculation')
    parser.add_argument('-n' ,'--nscf',        action="store_true", help='Non-self consistent calculation')
    parser.add_argument('-n2','--nscf_double', action="store_true", help='Non-self consistent calculation for the double grid')
    parser.add_argument('-b' ,'--bands',       action="store_true", help='Calculate band-structure')
    parser.add_argument('-p' ,'--phonon',      action="store_true", help='Phonon calculation')
    parser.add_argument('-t' ,'--nthreads',    action="store_true", help='Number of threads', default=2 )
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # create input files and folders
    relax()
    scf()
    nscf()
    bands()
   
    if args.relax: 
        print("running relax:")
        os.system("cd relax; mpirun -np %d pw.x -inp %s.scf > relax.log"%(args.nthreads,prefix))
        update_positions('relax','scf')
        print("done!")

    if args.scf:
        print("running scf:")
        os.system("cd scf; mpirun -np %d pw.x -inp %s.scf > scf.log"%(args.nthreads,prefix))
        print("done!")

    if args.nscf:
        print("running nscf:")
        os.system("cp -r scf/%s.save nscf/"%prefix)
        os.system("cd nscf; mpirun -np %d pw.x -inp %s.nscf > nscf.log"%(args.nthreads,prefix))
        print("done!")
    
    if args.bands:
        print("running bands:")
        os.system("cp -r scf/%s.save bands/"%prefix)
        os.system("cd bands; mpirun -np %d pw.x -inp %s.bands > bands.log"%(args.nthreads,prefix))
        print("done!")

        print("running plotting:")
        xml = PwXML(prefix='si',path='bands')
        xml.plot_eigen(p)
