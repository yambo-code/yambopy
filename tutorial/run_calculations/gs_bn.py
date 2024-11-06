#
# Author: Henrique Pereira Coutada Miranda
# Run a Silicon groundstate calculation using Quantum Espresso
#
from __future__ import print_function, division
import sys
import argparse
from qepy import *
from schedulerpy import *
from math import sqrt

kpoints      = [6,6,1]
kpoints_nscf = [6,6,1]
kpoints_double = [24,24,1]
qpoints = [3,3,1]
layer_separation = 12
pw = 'pw.x'
ph = 'ph.x'
q2r = 'q2r.x'
matdyn = 'matdyn.x'
prefix = 'bn'

npoints = 10 
p = Path([ [[0.0, 0.0, 0.0],r'$\Gamma$'],
           [[0.5, 0.0, 0.0],r'M'],
           [[1./3,1./3,0.0],r'K'],
           [[0.0, 0.0, 0.0],r'$\Gamma$']], [int(npoints*2),int(npoints),int(sqrt(5)*npoints)])

# scheduler
scheduler = Scheduler.factory

# create the input files
def get_inputfile():
    """ Define a Quantum espresso input file for boron nitride
    """ 
    qe = PwIn()
    qe.set_atoms([['N',[0.0,0.0,0.5]],
                  ['B',[1/3,2/3,0.5]]])
    qe.atypes = {'B': [10.811, "B.pbe-mt_fhi.UPF"],
                 'N': [14.0067,"N.pbe-mt_fhi.UPF"]}

    qe.control['prefix'] = "'%s'"%prefix
    qe.control['verbosity'] = "'high'"
    qe.control['wf_collect'] = '.true.'
    qe.control['pseudo_dir'] = "'../pseudos/'"
    qe.system['celldm(1)'] = 4.7
    qe.system['celldm(3)'] = layer_separation/qe.system['celldm(1)']
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 2
    qe.system['ibrav'] = 4
    qe.kpoints = [9, 9, 1]
    qe.electrons['conv_thr'] = 1e-10
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
    qe.write('relax/%s.relax'%prefix)

#scf
def scf(kpoints,folder='scf'):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    qe = get_inputfile()
    qe.control['calculation'] = "'scf'"
    qe.kpoints = kpoints
    qe.write('%s/%s.scf'%(folder,prefix))
 
#nscf
def nscf(kpoints,folder='nscf'):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    qe = get_inputfile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 70
    qe.system['force_symmorphic'] = ".true."
    qe.kpoints = kpoints
    qe.write('%s/%s.nscf'%(folder,prefix))

#bands
def bands():
    if not os.path.isdir('bands'):
        os.mkdir('bands')
    qe = get_inputfile()
    qe.control['calculation'] = "'bands'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-6
    qe.system['nbnd'] = 6
    qe.system['force_symmorphic'] = ".true."
    qe.ktype = 'crystal'
    qe.set_path(p)
    qe.write('bands/%s.bands'%prefix)

def phonon(kpoints,qpoints,folder='phonon'):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    ph = PhIn()
    ph['nq1'],ph['nq2'],ph['nq3'] = qpoints
    ph['tr2_ph'] = 1e-8
    ph['prefix'] = "'%s'"%prefix
    ph['epsil'] = ".false."
    ph['trans'] = ".true."
    ph['fildyn'] = "'%s.dyn'"%prefix
    ph['fildrho'] = "'%s.drho'"%prefix
    ph['ldisp'] = ".true."
    ph.write('%s/%s.ph'%(folder,prefix))

    md = DynmatIn()
    md['asr'] = "'simple'"
    md['fildyn'] = "'%s.dyn1'"%prefix
    md['filout'] = "'%s.modes'"%prefix
    md.write('%s/%s.dynmat'%(folder,prefix))

def update_positions(pathin,pathout):
    """ update the positions of the atoms in the scf file using the output of the relaxation loop
    """
    # Read scaled positions
    e = PwXML(prefix,path=pathin)
    pos = e.get_scaled_positions()
    
    #open relax input
    qin = PwIn.from_file('%s/%s.relax'%(pathin,prefix))
    print("old celldm(1)", qin.system['celldm(1)'])

    #open scf input
    qout = PwIn.from_file('%s/%s.scf'%(pathout,prefix))
    #replace lattice parameter
    qout.system['celldm(1)'] = e.cell[0][0]
    print("new celldm(1)", qout.system['celldm(1)'])
    #replace atomic positions
    new_atomic_pos = [[qout.atoms[i][0],list(pos[i])] for i in range(len(qout.atoms))]
    qout.set_atoms(new_atomic_pos)
    
    #re-write scf input
    qout.write('%s/%s.scf'%(pathout,prefix))

def run_plot():
    print("running plotting:")
    xml = PwXML(prefix=prefix,path='bands')
    xml.plot_eigen(p)

def run_projection(show=True):
    import matplotlib.pyplot as plt
    #write input file
    projwfc = ProjwfcIn('bn')
    projwfc.write(folder='bands')
    projwfc.run(folder='bands')
    #read xml file
    projection = ProjwfcXML(prefix='bn',path='bands')
    n_atom = range(16)
    b_atom = range(16,32)
    ax = plt.subplot(1,1,1)
    cax = projection.plot_eigen(ax,path=p,selected_orbitals=b_atom,selected_orbitals_2=n_atom,size=40,cmap='seismic')
    plt.colorbar(cax)
    if show: plt.show()

def run_bands(nthreads=1):
    print("running bands:")
    qe_run = scheduler() 
    qe_run.add_command("cp -r scf/%s.save bands/"%prefix)
    qe_run.add_command("cd bands; mpirun -np %d %s -inp %s.bands -nk %d > bands.log"%(nthreads,pw,prefix,nthreads))
    qe_run.run()
    qe_run.clean()
    print("done!")

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-r' ,'--relax',       action="store_true", help='Structural relaxation')
    parser.add_argument('-s' ,'--scf',         action="store_true", help='Self-consistent calculation')
    parser.add_argument('-n' ,'--nscf',        action="store_true", help='Non-self consistent calculation')
    parser.add_argument('-n2','--nscf_double', action="store_true", help='Non-self consistent calculation for the double grid')
    parser.add_argument('-b' ,'--bands',       action="store_true", help='Calculate band-structure')
    parser.add_argument('-l' ,'--plot',        action="store_true", help='Plot band-structure')
    parser.add_argument('-o' ,'--orbitals',    action="store_true", help='Plot atomic orbital projected band-structure')
    parser.add_argument('-p' ,'--phonon',      action="store_true", help='Phonon calculation')
    parser.add_argument('-d' ,'--dispersion',  action="store_true", help='Phonon dispersion')
    parser.add_argument('-t' ,'--nthreads',                         help='Number of threads', default=2 )
    args = parser.parse_args()
    nthreads = int(args.nthreads)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # create input files and folders
    relax()
    scf(kpoints,folder='scf')
    nscf(kpoints_nscf)
    nscf(kpoints_double, folder='nscf_double')
    bands()
    phonon(kpoints,qpoints)

    if args.relax:
        print("running relax:")
        qe_run = scheduler() 
        qe_run.add_command("cd relax; mpirun -np %d %s -inp %s.relax > relax.log"%(nthreads,pw,prefix))  #relax
        qe_run.run()
        update_positions('relax','scf') 
        print("done!")

    if args.scf:
        print("running scf:")
        qe_run = scheduler() 
        qe_run.add_command("cd scf; mpirun -np %d %s -inp %s.scf > scf.log"%(nthreads,pw,prefix))  #scf
        qe_run.run()
        print("done!")
   
    if args.nscf: 
        print("running nscf:")
        qe_run = scheduler() 
        qe_run.add_command("cp -r scf/%s.save nscf/"%prefix) #nscf
        qe_run.add_command("cd nscf; mpirun -np %d %s -nk %d -inp %s.nscf > nscf.log"%(nthreads,pw,nthreads,prefix)) #nscf
        qe_run.run()
        print("done!")

    if args.nscf_double: 
        print("running nscf_double:")
        qe_run = scheduler() 
        qe_run.add_command("cp -r scf/%s.save nscf_double/"%prefix) #nscf
        qe_run.add_command("cd nscf_double; mpirun -np %d %s -inp %s.nscf > nscf_double.log"%(nthreads,pw,prefix)) #nscf
        qe_run.run()
        print("done!")
    
    if args.phonon:
        print("running phonon:")
        qe_run = scheduler() 
        qe_run.add_command("cp -r scf/%s.save phonon/"%prefix)
        qe_run.add_command("cd phonon; mpirun -np %d %s -inp %s.ph > phonon.log"%(nthreads,ph,prefix)) #phonon
        qe_run.add_command("dynmat.x < %s.dynmat > dynmat.log"%prefix) #matdyn
        qe_run.run()
        print("done!")

    if args.dispersion:
        qe_run = scheduler() 

        #q2r
        disp = DynmatIn()
        disp['fildyn']= "'%s.dyn'" % prefix
        disp['zasr']  = "'simple'"
        disp['flfrc'] = "'%s.fc'"  % prefix
        disp.write('phonon/q2r.in')
        qe_run.add_command('cd phonon; %s < q2r.in'%q2r)

        #dynmat
        dyn = DynmatIn()
        dyn['flfrc'] = "'%s.fc'" % prefix
        dyn['asr']   = "'simple'"
        dyn['flfrq'] = "'%s.freq'" % prefix
        dyn['q_in_cryst_coord'] = '.true.'
        dyn.qpoints = p.get_klist()
        dyn.write('phonon/matdyn.in')
        qe_run.add_command('%s < matdyn.in'%matdyn)
        qe_run.run()

        # matdyn class to read and plot the frequencies
        m = Matdyn.from_modes_file(folder='phonon')
        m.plot_eigen(path=p)
 
    if args.bands:
        run_bands(nthreads)
        run_plot()

    if args.orbitals:
        run_projection()

    if args.plot:
        run_plot()

