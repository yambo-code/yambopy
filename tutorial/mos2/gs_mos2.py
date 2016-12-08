#
# Author: Henrique Pereira Coutada Miranda
# Run a MoS2 groundstate calculation using Quantum Espresso
#
from __future__ import print_function, division
from qepy import *
import argparse
import sys

scf_kpoints  = [12,12,1]
nscf_kpoints = [12,12,1]
dg_kpoints   = [24,24,1]
pw = 'pw.x'
ph = 'ph.x'
prefix = 'mos2'

npoints = 20
p = Path([ [[0.0, 0.0, 0.0],'G'],
           [[0.5, 0.0, 0.0],'M'],
           [[1./3,1./3,0.0],'K'],
           [[0.0, 0.0, 0.0],'G']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)])

#
# Create the input files
#
def get_inputfile():
    """ Define a Quantum espresso input file for MoS2
    """
    qe = PwIn()
    a = 5.838
    c = 20
    qe.atoms = [['Mo',[2/3,1/3,0.5]],
                [ 'S',[1/3,2/3, 2.92781466/c+0.5]],
                [ 'S',[1/3,2/3,-2.92781466/c+0.5]]]
    qe.atypes = {'Mo': [10.811, "Mo.pz-mt_fhi.UPF"],
                 'S':  [14.0067, "S.pz-mt_fhi.UPF"]}

    qe.control['prefix'] = "'mos2'"
    qe.control['wf_collect'] = '.true.'
    qe.control['verbosity'] = "'high'"
    qe.system['celldm(1)'] = a
    qe.system['celldm(3)'] = c/qe.system['celldm(1)']
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 3
    qe.system['ntyp'] = 2
    qe.system['ibrav'] = 4
    qe.kpoints = scf_kpoints
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
def nscf_kpoints_folder(kpoints,folder):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    qe = get_inputfile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 20
    qe.system['force_symmorphic'] = ".true."
    qe.kpoints = kpoints
    qe.write('%s/mos2.nscf'%folder)

#nscf
def nscf():
    nscf_kpoints_folder(nscf_kpoints,'nscf')

def nscf_double():
    nscf_kpoints_folder(dg_kpoints,'nscf_double')

#bands
def bands():
    if not os.path.isdir('bands'):
        os.mkdir('bands')
    qe = get_inputfile()
    qe.control['calculation'] = "'bands'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 13
    qe.system['force_symmorphic'] = ".true."
    qe.ktype = 'crystal'
    qe.set_path(p)
    qe.write('bands/%s.bands'%prefix)

def update_positions(pathin,pathout):
    """ update the positions of the atoms in the scf file using the output of the relaxation loop
    """
    e = PwXML('mos2',path=pathin)
    pos = e.get_scaled_positions()

    q = PwIn('%s/mos2.scf'%pathin)
    print("old celldm(1)", q.system['celldm(1)'])
    q.system['celldm(1)'] = e.cell[0][0]
    print("new celldm(1)", q.system['celldm(1)'])
    q.atoms = zip([a[0] for a in q.atoms],pos)
    q.write('%s/mos2.scf'%pathout)

def run_relax(nthreads=1):
    print("running relax:")
    os.system("cd relax; mpirun -np %d %s -inp mos2.scf > relax.log"%(nthreads,pw))  #relax
    update_positions('relax','scf')
    print("done!")

def run_scf(nthreads=1):
    print("running scf:")
    os.system("cd scf; mpirun -np %d %s -inp mos2.scf > scf.log"%(nthreads,pw))  #scf
    print("done!")

def run_nscf(nthreads=1):
    print("running nscf:")
    os.system("cp -r scf/mos2.save nscf/") #nscf
    os.system("cd nscf; mpirun -np %d %s -inp mos2.nscf -nk %d > nscf.log"%(nthreads,pw,nthreads)) #nscf
    print("done!")

def run_nscf_double(nthreads=1):
    print("running nscf_double:")
    os.system("cp -r scf/mos2.save nscf_double/") #nscf
    os.system("cd nscf_double; mpirun -np %d %s -inp mos2.nscf -nk %d > nscf_double.log"%(nthreads,pw,nthreads)) #nscf
    print("done!")

def run_bands(nthreads=1):
    print("running bands:")
    os.system("cp -r scf/%s.save bands/"%prefix)
    os.system("cd bands; mpirun -np %d %s -inp %s.bands -nk %d > bands.log"%(nthreads,pw,prefix,nthreads))
    print("done!")

def run_plot():
    print("running plotting:")
    xml = PwXML(prefix=prefix,path='bands')
    xml.plot_eigen(p)

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r' ,'--relax',       action="store_true", help='Structural relaxation')
    parser.add_argument('-s' ,'--scf',         action="store_true", help='Self-consistent calculation')
    parser.add_argument('-n' ,'--nscf',        action="store_true", help='Non self-consistent calculation')
    parser.add_argument('-n2','--nscf_double', action="store_true", help='Non self-consistent calculation for the double grid')
    parser.add_argument('-b' ,'--bands',       action="store_true", help='Calculate band-structure')
    parser.add_argument('-p' ,'--phonon',      action="store_true", help='Phonon calculation')
    parser.add_argument('-t' ,'--nthreads',    help='Number of threads', default=2 )
    args = parser.parse_args()

    nthreads = int(args.nthreads)
    print( "Using %d threads"%nthreads )

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # create input files and folders
    relax()
    scf()
    nscf()
    nscf_double()
    bands()

    if args.relax:
        run_relax(nthreads)

    if args.scf:
        run_scf(nthreads)

    if args.nscf:
        run_nscf(nthreads)

    if args.nscf_double:
        run_nscf_double(nthreads)

    if args.bands:
        run_bands(nthreads)
        run_plot()
