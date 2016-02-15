#
# Author: Henrique Pereira Coutada Miranda
# Run multiple GW calculations using Yambo for Boron Nitride
# using different layer separations with or without the coulomb truncation
#
from __future__ import division, print_function
from yambopy import *
from qepy import *
import argparse

yambo =  'yambo'

def get_inputfile(vac):
    """ Define a Quantum espresso input file for boron nitride
    """
    qe = PwIn()
    qe.atoms = [['B',[0.0,0.0,0.0]],
                ['N',[1/3,2/3,0.0]]]
    qe.atypes = {'B': [10.811, "B.pbe-mt_fhi.UPF"],
                 'N': [14.0067,"N.pbe-mt_fhi.UPF"]}

    qe.control['prefix'] = "'bn'"
    qe.control['wf_collect'] = '.true.'
    qe.system['celldm(1)'] = 4.7
    qe.system['celldm(3)'] = vac/qe.system['celldm(1)']
    qe.system['ecutwfc'] = 60
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 2
    qe.system['ibrav'] = 4
    qe.kpoints = [12, 12, 1]
    qe.electrons['conv_thr'] = 1e-8
    return qe

#scf
def scf(vac,folder):
    if not os.path.isdir('scf'):
        os.mkdir('scf')
    qe = get_inputfile(vac)
    qe.control['calculation'] = "'scf'"
    qe.write('%s/bn.scf'%folder)

#nscf
def nscf(vac,kpoints,folder):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    qe = get_inputfile(vac)
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = 30
    qe.system['force_symmorphic'] = ".true."
    qe.kpoints = kpoints
    qe.write('%s/bn.nscf'%folder)

#parse options
parser = argparse.ArgumentParser(description='Test the yambopy script.')
parser.add_argument('-r','--run',      action="store_true", help='Run structural relaxation')
parser.add_argument('-a','--analyse',  action="store_true", help='Run non-self consistent calculation for the double grid')
parser.add_argument('-c','--cut',      action="store_true", help='Run non-self consistent calculation for the double grid')
parser.add_argument('-t','--nthreads', action="store_true", help='Number of threads', default=2 )
args = parser.parse_args()

separations = [12,14,16,18,20,22,23,24,25]
if args.run:
    for vac in separations:

        folder = 'gw_cutoff/%d'%vac
        scf_folder  = '%s/scf'%folder
        nscf_folder = '%s/nscf'%folder
        os.system('mkdir -p %s'%scf_folder)
        os.system('mkdir -p %s'%nscf_folder)

        #if database not present calculate it
        if not os.path.isdir("%s/SAVE"%folder):
            print("vaccum: %d"%vac)
            print('calculate scf')
            scf(vac,'%s/scf'%folder)
            os.system("cd %s; mpirun -np %d pw.x -inp bn.scf > scf.log"%(scf_folder,args.nthreads))  #scf

            print('calculate nscf')
            os.system('cp -r %s/bn.save %s/'%(scf_folder,nscf_folder))
            nscf(vac,[6,6,1],'%s/nscf'%folder)
            os.system("cd %s; mpirun -np %d pw.x -inp bn.nscf > nscf.log"%(nscf_folder,args.nthreads)) #nscf

            print('run p2y and yambo')
            os.system('cd %s/bn.save; p2y > p2y.log'%nscf_folder)
            os.system('cd %s/bn.save; yambo > yambo.log'%nscf_folder)
            os.system('mv %s/bn.save/SAVE %s'%(nscf_folder,folder))

        #create the yambo input file
        if args.cut:
            y = YamboIn('%s -r -d -g n -V all'%yambo,folder=folder)
        else:
            y = YamboIn('%s -d -g n -V all'%yambo,folder=folder)

        QPKrange,_ = y['QPkrange']
        y['QPkrange'] = [QPKrange[:2]+[3,6],'']
        y['FFTGvecs'] = [15,'Ry']
        y['NGsBlkXd'] = [1,'RL']
        y['BndsRnXd'] = [[1,30],'']
        y['ETStpsXd'] = [50,'']

        if args.cut:
            y['CUTGeo'] = "box z"                    # [CUT] Coulomb Cutoff geometry: box/cylinder/sphere X/Y/Z/XY..
            y['CUTBox'] = [[ 0.00, 0.00, vac-2],'' ] # [CUT] [au] Box sides

        y.write('gw_cutoff/%d/yambo_run.in'%vac)

        print('running yambo')
        os.system('cd %s; %s -F yambo_run.in -J vac_%d -C yambo > vac_%d.log'%(folder,yambo,vac,vac))

if args.analyse:
    #collect all the data pack the files in .json files
    for vac in separations:
        folder = 'gw_cutoff/%d/yambo'%vac
        y = YamboOut(folder)
        if not y.locked():
            y.pack()
            y.put_lock()

        os.system('cp gw_cutoff/%d/yambo.json gw_cutoff/%d.json'%(vac,vac))

    #plot the band structure
    ya = YamboAnalyser('gw_cutoff')
    ya.plot_gw(['qp'],cols=(lambda x: x[2]+x[3],))
    ya.plot_gw(['qp'],cols=(lambda x: x[2]+x[3],),rows=(lambda x: x[2]-x[1],))
