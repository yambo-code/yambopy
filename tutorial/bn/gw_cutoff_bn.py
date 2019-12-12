#
# Author: Henrique Pereira Coutada Miranda
# Run multiple GW calculations using Yambo for Boron Nitride
# using different layer separations with or without the coulomb truncation
#
from __future__ import division, print_function
from yambopy import *
from qepy import *
from schedulerpy import *
import argparse

yambo =  'yambo'

# scheduler
scheduler = Scheduler.factory

def get_inputfile(vac):
    """ Define a Quantum espresso input file for boron nitride
    """
    qe = PwIn()
    qe.set_atoms([['N',[0.0,0.0,0.5]],
                  ['B',[1/3,2/3,0.5]]])
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
parser.add_argument('-r','--run',      action="store_true", help='Run scf+nscf+GW calculations for different vacuum distances')
parser.add_argument('-a','--analyse',  action="store_true", help='Analyse GW results')
parser.add_argument('-c','--cut',      action="store_true", help='Use Coulomb cutoff in GW runs')
parser.add_argument('-t','--nthreads', action="store_true", help='Number of threads', default=2 )
args = parser.parse_args()

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

separations = [12,14,16,18,20,22,23,24,25]
if args.run:
    for vac in separations:
        shell=scheduler()
        folder = 'gw_cutoff/%d'%vac
        scf_folder  = '%s/scf'%folder
        nscf_folder = '%s/nscf'%folder
        shell.add_command('mkdir -p %s'%scf_folder)
        shell.add_command('mkdir -p %s'%nscf_folder)
        shell.run()

        #if database not present calculate it
        if not os.path.isdir("%s/SAVE"%folder):
            scf_run = scheduler()
            print("vacuum: %d"%vac)
            print('calculate scf')
            scf(vac,'%s/scf'%folder)
            scf_run.add_command("cd %s; mpirun -np %d pw.x -inp bn.scf > scf.log"%(scf_folder,args.nthreads))  #scf
            scf_run.run()

            print('calculate nscf')
            nscf_run = scheduler()
            nscf_run.add_command('cp -r %s/bn.save %s/'%(scf_folder,nscf_folder))
            nscf(vac,[6,6,1],'%s/nscf'%folder)
            nscf_run.add_command("cd %s; mpirun -np %d pw.x -inp bn.nscf > nscf.log"%(nscf_folder,args.nthreads)) #nscf
            nscf_run.run()

            print('run p2y and yambo')
            p2y_run = scheduler()
            p2y_run.add_command('cd %s/bn.save; p2y > p2y.log'%nscf_folder)
            p2y_run.add_command('yambo > yambo.log')
            p2y_run.add_command('mv SAVE ../../')
            p2y_run.run()

        #create the yambo input file
        if args.cut:
            y = YamboIn.from_runlevel('%s -r -d -g n -V all'%yambo,folder=folder)
        else:
            y = YamboIn.from_runlevel('%s -d -g n -V all'%yambo,folder=folder)

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
        run_yambo = scheduler()
        run_yambo.add_command('cd %s; %s -F yambo_run.in -J vac_%d -C yambo > vac_%d.log'%(folder,yambo,vac,vac))
        run_yambo.run()

if args.analyse:
    #collect all the data pack the files in .json files
    for vac in separations:
        folder = 'gw_cutoff/%d/yambo'%vac
        save_folder = 'gw_cutoff/%d'%vac
        y = YamboOut(folder=folder,save_folder=save_folder)
        y.pack()    

        shell = scheduler()
        shell.add_command('cp gw_cutoff/%d/yambo.json gw_cutoff/%d.json'%(vac,vac))
        shell.run()

    #plot the band structure
    ya = YamboAnalyser('gw_cutoff')
    ya.plot_gw(['qp'],cols=(lambda x: x[2]+x[3],))
    ya.plot_gw(['qp'],cols=(lambda x: x[2]+x[3],),rows=(lambda x: x[2]-x[1],))
