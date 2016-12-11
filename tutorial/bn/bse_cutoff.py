#
# Author: Henrique Pereira Coutada Miranda
# Check the convergence of the coulomb cutoff for a BSE calculation using yambo
#
from __future__ import print_function
from yambopy import *
from qepy import *
import argparse
import sys

prefix = "bn"
work_folder = "bse_cutoff"
yambo = "yambo"
layer_separations = [10,12,14,16]
kpoints = [9,9,1]

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
    qe.system['celldm(3)'] = 12/qe.system['celldm(1)']
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 2
    qe.system['ibrav'] = 4
    qe.kpoints = [9, 9, 1]
    qe.electrons['conv_thr'] = 1e-8
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
    qe.system['nbnd'] = 30
    qe.system['force_symmorphic'] = ".true."
    qe.system['celldm(3)'] = layer_separation/qe.system['celldm(1)']
    qe.kpoints = kpoints
    qe.write('%s/%s.nscf'%(folder,prefix))
    
def database(output_folder,nscf_folder='nscf'):
    if not os.path.isdir('%s/SAVE'%output_folder):
        print('preparing yambo database...')
        os.system('cd %s/%s.save; p2y > p2y.log'%(nscf_folder,prefix))
        os.system('cd %s/%s.save; yambo > yambo.log'%(nscf_folder,prefix))
        os.system('mv %s/%s.save/SAVE %s'%(nscf_folder,prefix,output_folder))
        print('done!')

def run(nthreads=1):
    #for each separation run the ground state calculation and
    for layer_separation in layer_separations:

        print("layer separation: %d bohr"%layer_separation)
        root_folder = "%s/%d"%(work_folder,layer_separation)
        if not os.path.isdir(root_folder):
            os.makedirs(root_folder)

        # run the ground state calculation
        print("scf cycle")
        scf(layer_separation,folder="%s/scf"%root_folder)
        os.system("cd %s/scf; pw.x < %s.scf > scf.log"%(root_folder,prefix))

        # run the non self consistent calculation
        print("nscf cycle")
        src ='%s/scf/%s.save'%(root_folder,prefix)
        dst ='%s/nscf/%s.save'%(root_folder,prefix)
        nscf(layer_separation,folder="%s/nscf"%root_folder)
        os.system( 'cp -r %s %s'%(src,dst) )
        os.system("cd %s/nscf; pw.x < %s.nscf > nscf.log"%(root_folder,prefix))

        # generate the database
        database('%s'%root_folder,nscf_folder="%s/nscf"%root_folder)

        # calculate the absorption spectra
        y = YamboIn('mpirun -np %d yambo -r -b -o b -k sex -y d -V all'%nthreads,folder=root_folder)

        y['FFTGvecs'] = [30,'Ry']
        y['NGsBlkXs'] = [1,'Ry']
        y['BndsRnXs'] = [1,30]

        y['CUTGeo'] = 'box z'
        y['CUTBox'] = [0,0,layer_separation-1]

        y['KfnQP_E']  = [1.0,1.0,1.0] #scissor operator
        y['BSEBands'] = [3,6]
        y['BEnSteps'] = 500
        y['BEnRange'] = [[1.0,6.0],'eV']
        y.write('%s/yambo_run.in'%root_folder)
        os.system('cd %s; %s -F yambo_run.in -J %d'%(root_folder,yambo,layer_separation))

def plot():
    for layer_separation in layer_separations:
        root_folder = "%s/%d"%(work_folder,layer_separation)
    
        #gather the results
        pack_files_in_folder(root_folder)

    #plot the results
    ya = YamboAnalyser(work_folder)
    ya.plot_bse('eps')

if __name__ == "__main__":

    #parse options
    parser = argparse.ArgumentParser(description='Convergence test of the colomb cutoff')
    parser.add_argument('-r' ,'--run',     action="store_true", help='Run the calculation')
    parser.add_argument('-p' ,'--plot',    action="store_true", help='Run the analysis')
    parser.add_argument('-t' ,'--nthreads',                     help='Number of threads', default=2)
    args = parser.parse_args()
    nthreads = int(args.nthreads)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    if args.run:
        run(nthreads)
    if args.plot:
        plot()
