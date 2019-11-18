#
# Author: Henrique Pereira Coutada Miranda
# Run a Silicon groundstate calculation using Quantum Espresso
#
from __future__ import print_function
import os
import sys
import argparse
from qepy import *

scf_kpoints  = [2,2,2]
nscf_kpoints = [2,2,2]
dg_kpoints   = [4,4,4]
prefix = 'si'
matdyn = 'matdyn.x'
q2r =    'q2r.x'
pw = 'pw.x'
ph = 'ph.x'
p = Path([ [[1.0,1.0,1.0],'$\Gamma$'],
           [[0.0,0.5,0.5],'$X$'],
           [[0.0,0.0,0.0],'$\Gamma$'],
           [[0.5,0.0,0.0],'$L$']], [20,20,20])

#
# Create the input files
#
def get_inputfile():
    """ Define a Quantum espresso input file for silicon
    """
    qe = PwIn()
    qe.set_atoms([['Si',[0.125,0.125,0.125]],
                  ['Si',[-.125,-.125,-.125]]])
    qe.atypes = {'Si': [28.086,"Si.pbe-mt_fhi.UPF"]}

    qe.control['prefix'] = "'%s'"%prefix
    qe.control['wf_collect'] = '.true.'
    qe.control['pseudo_dir'] = "'../pseudos'"
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
    qe.system['nbnd'] = 30
    qe.system['force_symmorphic'] = ".true."
    qe.kpoints = nscf_kpoints
    qe.write('nscf/%s.nscf'%prefix)

#double-grid
def dg():
    if not os.path.isdir('nscf-dg'):
        os.mkdir('nscf-dg')
    qe = get_inputfile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 8
    qe.system['force_symmorphic'] = ".true."
    qe.kpoints = dg_kpoints
    qe.write('nscf-dg/%s.nscf'%prefix)


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

def plot_orbitals(show=True):
    import matplotlib.pyplot as plt
    projwfc = ProjwfcIn('si')
    projwfc.write(folder='bands')
    projwfc.run(folder='bands')
    projection = ProjwfcXML(prefix='si',path='bands')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    s_orb = [0,16]
    p_orb = [1,2,3,17,19,20]
    projection.plot_eigen(ax,path=p,selected_orbitals=s_orb,selected_orbitals_2=p_orb,size=40,cmap='RdBu')
    ax.set_ylim([-7,6])
    if show: plt.show()

def phonons():
    os.system('mkdir -p phonons')
    ph = PhIn()
    ph['prefix'] = "'%s'"     % prefix
    ph['fildyn'] = "'%s.dyn'" % prefix
    ph['ldisp']  = '.true.'
    ph['trans']  = '.true.'
    ph['tr2_ph'] = 1e-12
    ph['nq1'], ph['nq2'], ph['nq3'] = 2, 2, 2
    ph.write('phonons/%s.phonons'%prefix)

def dispersion():
    disp = DynmatIn()
    disp['fildyn']= "'%s.dyn'" % prefix
    disp['zasr']  = "'simple'" 
    disp['flfrc'] = "'%s.fc'"  % prefix
    disp.write('phonons/q2r.in')
    os.system('cd phonons; %s < q2r.in'%q2r)
    dyn = DynmatIn()
    dyn['flfrc'] = "'%s.fc'" % prefix
    dyn['asr']   = "'simple'"  
    dyn['flfrq'] = "'%s.freq'" % prefix
    dyn['q_in_cryst_coord'] = '.true.'
    dyn.qpoints = p.get_klist()
    dyn.write('phonons/matdyn.in')
    os.system('cd phonons; %s < matdyn.in'%matdyn)
    print( len(p.get_klist()) ) 
    # Use a class to read and plot the frequencies
    Matdyn(natoms=2,path=p,folder='phonons').plot_eigen()

def update_positions(pathin,pathout):
    """ update the positions of the atoms in the scf file using the output of the relaxation loop
    """
    e = PwXML(prefix,path=pathin)
    pos = e.get_scaled_positions()
     
    #open relaxed cell
    qin  = PwIn('%s/%s.scf'%(pathin,prefix))
  
    #open scf file
    qout = PwIn('%s/%s.scf'%(pathout,prefix))

    #update positions on scf file
    print("old celldm(1)", qin.system['celldm(1)'])
    qout.system['celldm(1)'] = e.cell[0][2]*2
    print("new celldm(1)", qout.system['celldm(1)'])
    qout.atoms = list(zip([a[0] for a in qin.atoms],pos))

    #write scf
    qout.write('%s/%s.scf'%(pathout,prefix))

def run_relax(nthreads=1):
    print("running relax:")
    os.system("cd relax; mpirun -np %d %s -inp %s.scf > relax.log"%(nthreads,pw,prefix))
    update_positions('relax', 'scf')
    print("done!")

def run_scf(nthreads=1):
    print("running scf:")
    os.system("cd scf; mpirun -np %d %s -inp %s.scf > scf.log"%(nthreads,pw,prefix))
    print("done!")

def run_nscf(nthreads=1):
    print("running nscf:")
    os.system("cp -r scf/%s.save nscf/"%prefix)
    os.system("cd nscf; mpirun -np %d %s -inp %s.nscf > nscf.log"%(nthreads,pw,prefix))
    print("done!")

def run_dg(nthreads=1):
    print("running nscf:")
    os.system("cp -r scf/%s.save nscf-dg/"%prefix)
    os.system("cd nscf-dg; mpirun -np %d %s -inp %s.nscf > nscf.log"%(nthreads,pw,prefix))
    print("done!")

def run_bands(nthreads=1):
    print("running bands:")
    os.system("cp -r scf/%s.save bands/"%prefix)
    os.system("cd bands; mpirun -np %d %s -inp %s.bands > bands.log"%(nthreads,pw,prefix))
    print("done!")

def run_plot(show=True):
    print("running plotting:")
    xml = PwXML(prefix='si',path='bands')
    xml.plot_eigen(p,show=show)

def run_phonon(nthreads=1):
    print("running phonons:")
    os.system("cp -r scf/%s.save phonons/"%prefix)
    os.system("cd phonons; mpirun -np %d %s -inp %s.phonons > phonons.log"%(nthreads,ph,prefix))
    print("done!")

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r' ,'--relax',       action="store_true", help='Structural relaxation')
    parser.add_argument('-s' ,'--scf',         action="store_true", help='Self-consistent calculation')
    parser.add_argument('-n' ,'--nscf',        action="store_true", help='Non-self consistent calculation')
    parser.add_argument('-n2','--nscf_double', action="store_true", help='Non-self consistent calculation for the double grid')
    parser.add_argument('-b' ,'--bands',       action="store_true", help='Calculate band-structure')
    parser.add_argument('-o' ,'--orbitals',    action="store_true", help='Plot band structure with orbital weights')
    parser.add_argument('-p' ,'--phonon',      action="store_true", help='Phonon calculation')
    parser.add_argument('-d' ,'--dispersion',  action="store_true", help='Phonon dispersion')
    parser.add_argument('-t' ,'--nthreads',                         help='Number of threads', default=2 )
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # create input files and folders
   
    scf()
    if args.relax:
        relax()
        run_relax(args.nthreads) 
    if args.scf:        
        run_scf(args.nthreads)
    if args.nscf:
        nscf()
        run_nscf(args.nthreads)
    if args.nscf_double:
        dg()
        run_dg(args.nthreads)
    if args.phonon:     
        phonons()
        run_phonon(args.nthreads)
    if args.dispersion: dispersion()
    if args.bands:
        bands()
        run_bands(args.nthreads)
        run_plot()
    if args.orbitals:
        plot_orbitals()

