#
# Author: Henrique Pereira Coutada Miranda
# Run a Silicon groundstate calculation using Quantum Espresso
#
from __future__ import print_function
from pwpy.inputfile import *
from pwpy.outputxml import *

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

    qe.control['prefix'] = "'si'"
    qe.control['wf_collect'] = '.true.'
    qe.system['celldm(1)'] = 10.3
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 1
    qe.system['ibrav'] = 2
    qe.kpoints = [4, 4, 4]
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
    qe.write('relax/si.scf')

#scf
def scf():
    if not os.path.isdir('scf'):
        os.mkdir('scf')
    qe = get_inputfile()
    qe.control['calculation'] = "'scf'"
    qe.write('scf/si.scf')
 
#nscf
def nscf():
    if not os.path.isdir('nscf'):
        os.mkdir('nscf')
    qe = get_inputfile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 30
    qe.system['force_symmorphic'] = ".true."
    qe.kpoints = [2, 2, 2]
    qe.write('nscf/si.nscf')

def update_positions(pathin,pathout):
    """ update the positions of the atoms in the scf file using the output of the relaxation loop
    """
    e = EspressoXML('si',path=pathin)
    pos = e.get_scaled_positions()

    q = PwIn('%s/si.scf'%pathin)
    print("old celldm(1)", q.system['celldm(1)'])
    q.system['celldm(1)'] = e.cell[0][2]*2
    print("new celldm(1)", q.system['celldm(1)'])
    q.atoms = zip([a[0] for a in q.atoms],pos)
    q.write('%s/si.scf'%pathout)

if __name__ == "__main__":
    nproc = 1

    # create input files and folders
    relax()
    scf()
    nscf()

    print("running relax:")
    os.system("cd relax; mpirun -np %d pw.x -inp si.scf > relax.log"%nproc)
    update_positions('relax','scf') 
    print("done!")

    print("running scf:")
    os.system("cd scf; mpirun -np %d pw.x -inp si.scf > scf.log"%nproc)
    print("done!")

    print("running nscf:")
    os.system("cp -r scf/si.save nscf/")
    os.system("cd nscf; mpirun -np %d pw.x -inp si.nscf > nscf.log"%nproc)
    print("done!")

