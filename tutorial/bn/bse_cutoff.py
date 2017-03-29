#
# Author: Henrique Pereira Coutada Miranda
# Check the convergence of the coulomb cutoff for a BSE calculation using yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *
from schedulerpy import *
from functools import partial
import multiprocessing
import argparse
import sys

prefix = "bn"
yambo = "yambo"
p2y = 'p2y'
layer_separations = [10,15,20,25,30,35,40] #,45,50]
kpoints = [24,24,1]
nbands = 20 #should be larger! we use 30 for demonstration

scheduler = Scheduler.factory

# create the quantum espresso input file
def get_inputfile():
    """ Define a Quantum espresso input file for boron nitride
    """ 
    qe = PwIn()
    qe.atoms = [['N',[ 0.0, 0.0,0.5]],
                ['B',[1./3,2./3,0.5]]]
    qe.atypes = {'B': [10.811, "B.pbe-mt_fhi.UPF"],
                 'N': [14.0067,"N.pbe-mt_fhi.UPF"]}

    qe.control['prefix'] = "'%s'"%prefix
    qe.control['wf_collect'] = '.true.'
    qe.system['celldm(1)'] = 4.7
    qe.system['celldm(3)'] = 14/qe.system['celldm(1)']
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 2
    qe.system['ibrav'] = 4
    qe.kpoints = [9, 9, 1]
    qe.electrons['conv_thr'] = 1e-10
    return qe

#run the self consistent calculation
def scf(layer_separation,folder='scf'):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    qe = get_inputfile()
    qe.system['celldm(3)'] = layer_separation/qe.system['celldm(1)']
    qe.control['calculation'] = "'scf'"
    qe.write('%s/%s.scf'%(folder,prefix))

#run the non-self consistent calculation
def nscf(layer_separation,folder='nscf'):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    qe = get_inputfile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = nbands
    qe.system['force_symmorphic'] = ".true."
    qe.system['celldm(3)'] = layer_separation/qe.system['celldm(1)']
    qe.kpoints = kpoints
    qe.write('%s/%s.nscf'%(folder,prefix))
    
def database(shell,output_folder,nscf_folder='nscf'):
    if not os.path.isdir('%s/SAVE'%output_folder):
        print('preparing yambo database...')
        shell.add_command('mkdir -p %s'%nscf_folder)
        shell.add_command('pushd %s/%s.save; %s; %s'%(nscf_folder,prefix,p2y,yambo))
        shell.add_command('popd')
        shell.add_command('mv %s/%s.save/SAVE %s'%(nscf_folder,prefix,output_folder))
        print('done!')

def run_job(layer_separation,work_folder='bse_cutoff',cut=False):
    """
    Given a layer separation run the calculation
    1. scf calculation with QE
    2. nscf calculation
    3. BSE with yambo
    """
    #check if the calculation exists
    done_stamp = '%s/%d/done'%(work_folder,layer_separation)
    print(done_stamp)
    if os.path.isfile(done_stamp):
        return

    print("layer separation: %d bohr    cutoff:"%layer_separation, cut)
    root_folder = "%s/%d"%(work_folder,layer_separation)
    shell = scheduler()
    if not os.path.isdir(root_folder):
        shell.add_command( 'mkdir -p %s'%root_folder )

    # 1. run the ground state calculation
    print("scf cycle")
    scf(layer_separation,folder="%s/scf"%root_folder)
    shell.add_command("pushd %s/scf; pw.x < %s.scf > scf.log"%(root_folder,prefix))
    shell.add_command("popd")

    # 2. run the non self consistent calculation
    print("nscf cycle")
    src ='%s/scf/%s.save'%(root_folder,prefix)
    dst ='%s/nscf/%s.save'%(root_folder,prefix)
    nscf(layer_separation,folder="%s/nscf"%root_folder)

    shell.add_command('cp -r %s %s'%(src,dst) )
    shell.add_command("pushd %s/nscf; pw.x < %s.nscf > nscf.log"%(root_folder,prefix))
    shell.add_command('popd')

    # generate the database
    database(shell,'%s'%root_folder,nscf_folder="%s/nscf"%root_folder)
    shell.run()
    #wait for execution

    # 3. calculate the absorption spectra
    y = YamboIn('mpirun -np %d yambo -r -b -o b -k sex -y d -V all'%nthreads,folder=root_folder)

    if cut:
        y['CUTGeo'] = 'box z'
        y['CUTBox'] = [0,0,layer_separation-2]

        y['RandQpts'] = 1000000
        y['RandGvec'] = [1,'Ry']

    y['FFTGvecs'] = [20,'Ry']
    y['NGsBlkXs'] = [1,'Ry']          #local field effects
    y['BndsRnXs'] = [1,nbands]            #number of bands for static screening

    y['KfnQP_E']  = [2.91355133,1.0,1.0] #scissor operator
    y['BSEBands'] = [4,5]                #number of bands in BSE kernel
    y['BEnRange'] = [[4.0,8.0],'eV']     #energy range to plot optical absorption
    y['BEnSteps'] = 500                  #energy steps in the range
    y.write('%s/yambo_run.in'%root_folder)

    shell = scheduler()
    shell.add_command('cd %s; %s -F yambo_run.in -J %d'%(root_folder,yambo,layer_separation))
    shell.add_command('touch done')
    shell.run()

def run(nthreads=1,work_folder='bse_cutoff',cut=True):

    if (nthreads > 1):
        p = multiprocessing.Pool(nthreads)
        run = partial(run_job,work_folder=work_folder,cut=cut)
        try:
            #reversed list because of load imbalance
            p.map(run, reversed(layer_separations))
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            p.terminate()
            p.join()

    else:
        for layer_separation in layer_separations:
            run_job(layer_separation,work_folder,cut)

def plot(work_folder):
    ax = plt.gca()
    for layer_separation in layer_separations:
        root_folder = "%s/%d"%(work_folder,layer_separation)
    
        #gather the results
        pack_files_in_folder(root_folder)

    #plot the results
    ya = YamboAnalyser(work_folder)
    print(ya)
    ax = ya.plot_bse('eps',ax=ax)

    if cut: title = "with coulomb cutoff"
    else:   title = "without coulomb cutoff"

    plt.title(title)
    plt.show()

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Convergence test of the colomb cutoff')
    parser.add_argument('-r' ,'--run',     action="store_true", help='Run the calculation')
    parser.add_argument('-c' ,'--cut',     action="store_true", help='Use coulomb cutoff')
    parser.add_argument('-p' ,'--plot',    action="store_true", help='Run the analysis')
    parser.add_argument('-t' ,'--nthreads',                     help='Number of threads', default=1)
    args = parser.parse_args()
    nthreads = int(args.nthreads)
    print("Using %d threads"%nthreads)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    cut = args.cut
    
    #choose work_folder
    if cut: 
        work_folder = "bse_cutoff_cut"
    else:
        work_folder = "bse_cutoff"

    if args.run:
        run(nthreads,work_folder,cut)
    if args.plot:
        plot(work_folder)
