#
#
# Tutorial Yambo School. Lausanne, 24-28 April 2017
# Convergence GW on hexagonal BN
# Alejandro Molina-Sanchez & Henrique P. C. Miranda
#
#
from __future__ import print_function
import sys
from yambopy     import *
from qepy        import *
from schedulerpy import *
import argparse

yambo = 'yambo'
p2y = 'p2y'
prefix = 'bn'

bash = Scheduler.factory

def create_save():
    #check if the nscf cycle is present
    if os.path.isdir('nscf/%s.save'%prefix):
        print('nscf calculation found!')
    else:
        print('nscf calculation not found!')
        exit() 

    #check if the SAVE folder is present
    if not os.path.isdir('database'):
        print('preparing yambo database')
        shell = bash() 
        shell.add_command('mkdir -p database')
        shell.add_command('cd nscf/%s.save; %s; %s'%(prefix, p2y, yambo))
        shell.add_command('mv SAVE  ../../database/')
        shell.run()
        shell.clean()

def gw_convergence():
    #create the folder to run the calculation
    if not os.path.isdir('gw_conv'):
        shell = bash() 
        shell.add_command('mkdir -p gw_conv')
        shell.add_command('cp -r database/SAVE gw_conv/')
        shell.run()
        shell.clean()

    y = YamboIn('%s -p p -g n -V all'%yambo,folder='gw_conv')
    k_f = y['QPkrange'][0][1]         # Read the last k-points in the uniform k-grid

    y['BndsRnXp'] = [[1,10],'']             # Screening. Number of bands
    y['NGsBlkXp'] = [0,'Ry']                # Cutoff Screening
    y['GbndRnge'] = [[1,10],'']             # Self-energy. Number of bands
    y['QPkrange'] = [ [k_f,k_f,4,5], '' ]

    conv = { 'EXXRLvcs': [[10,10,20,40,60,80,100,150],'Ry'],
             'NGsBlkXp': [[0,0,1,2,3], 'Ry'],
             'BndsRnXp': [[[1,10],[1,10],[1,15],[1,20],[1,30]],''] ,
             'GbndRnge': [[[1,10],[1,10],[1,15],[1,20],[1,30]],''] }

    def run(filename):
        """ Function to be called by the optimize function """
        folder = filename.split('.')[0]
        print(filename,folder)
        shell = bash() 
        shell.add_command('cd gw_conv')
        shell.add_command('rm -f *.json %s/o-*'%folder) #cleanup
        shell.add_command('%s -F %s -J %s -C %s 2> %s.log'%(yambo,filename,folder,folder,folder))
        shell.run()
        shell.clean()

    y.optimize(conv,run=run,ref_run=False)

def plot_convergence():
    y = YamboIn('%s -d -g n -V all'%yambo,folder='gw_conv')

    k_f = y['QPkrange'][0][1]         # Read last k-points in the uniform k-grid
    print (k_f) #pack the files in .json files
    pack_files_in_folder('gw_conv')

    print('Select the converged value for each variable')
    shell = bash() 
    shell.add_command('yambopy analysegw -bc 5 -kc %s -bv 4 -kv %s gw_conv EXXRLvcs' % (k_f, k_f))
    shell.add_command('yambopy analysegw -bc 5 -kc %s -bv 4 -kv %s gw_conv NGsBlkXp' % (k_f, k_f))
    shell.add_command('yambopy analysegw -bc 5 -kc %s -bv 4 -kv %s gw_conv BndsRnXp' % (k_f, k_f))
    shell.add_command('yambopy analysegw -bc 5 -kc %s -bv 4 -kv %s gw_conv GbndRnge' % (k_f, k_f))
    shell.run()
    shell.clean()

def gw():
    #create the folder to run the calculation
    if not os.path.isdir('gw'):
        shell = bash() 
        shell.add_command('mkdir -p gw')
        shell.add_command('cp -r database/SAVE gw/')
        shell.run()
        shell.clean()

    # GW calculation. PPA Screening. Newton method
    y = YamboIn('%s -p p -g n -V all'%yambo,folder='gw')

    y['EXXRLvcs'] = [80,'Ry']       # Self-energy. Exchange
    y['BndsRnXp'] = [1,25]          # Screening. Number of bands
    y['NGsBlkXp'] = [3,'Ry']        # Cutoff Screening
    y['GbndRnge'] = [1,25]          # Self-energy. Number of bands
    #read values from QPkrange
    values, units = y['QPkrange']
    kpoint_start, kpoint_end, band_start, band_end = values
    #set the values of QPkrange
    y['QPkrange'] = [kpoint_start,kpoint_end,2,6]
    y.write('gw/yambo_gw.in')

    print('calculating...')
    shell = bash() 
    shell.add_command('cd gw')
    shell.add_command('rm -f *.json gw/o-*') #cleanup
    shell.add_command('%s -F yambo_gw.in -J gw -C gw' % yambo)
    shell.run()
    shell.clean()

def plot_gw():
    #pack the files in .json files
    pack_files_in_folder('gw')

    #plot the results using yambm analyser
    ya = YamboAnalyser('gw')
    print('plot all qpoints')
    ya.plot_gw('qp')

    print('plot along a path')
    path = [[[0,   0,   0],'$\Gamma$'],
            [[0.5, 0,   0],'M'],
            [[0.3333,0.3333, 0.0],'K'],
            [[0.0, 0.0, 0.0],'$\Gamma$']]
    ya.plot_gw_path('qp',path, cols=(lambda x: x[2]+x[3],2))

def xi():
    #create the folder to run the calculation
    if not os.path.isdir('gw-xi'):
        shell = bash() 
        shell.add_command('mkdir -p gw-xi')
        shell.add_command('cp -r database/SAVE gw-xi/')
        shell.run()
        shell.clean()

    print ("Running COHSEX in folder 'gw-xi/coh'")
    cohsex = YamboIn('%s -p c -g n -V all'%yambo,folder='gw-xi')
    cohsex['EXXRLvcs'] = [80,'Ry']       # Self-energy. Exchange
    cohsex['BndsRnXs'] = [1,25]          # Screening. Number of bands
    cohsex['NGsBlkXs'] = [3,'Ry']        # Cutoff Screening
    cohsex['GbndRnge'] = [1,25]          # Self-energy. Number of bands
    cohsex['QPkrange'][0][2:] = [2,6]
    cohsex.write('gw-xi/yambo_cohsex.in')
    shell = bash() 
    shell.add_command('cd gw-xi')
    shell.add_command('rm -f coh.json coh/o-coh*') #cleanup
    shell.add_command('%s -F yambo_cohsex.in -J coh -C coh' % yambo)
    shell.run()
    shell.clean()

    print ("Running COHSEX in folder 'gw-xi/pp'")
    ppa = YamboIn('%s -p p -g n -V all'%yambo,folder='gw-xi')
    ppa['EXXRLvcs'] = [80,'Ry']       # Self-energy. Exchange
    ppa['BndsRnXp'] = [1,25]          # Screening. Number of bands
    ppa['NGsBlkXp'] = [3,'Ry']        # Cutoff Screening
    ppa['GbndRnge'] = [1,25]          # Self-energy. Number of bands
    ppa['QPkrange'][0][2:] = [2, 6]       # QP range. All BZ
    ppa.write('gw-xi/yambo_ppa.in')
    shell = bash() 
    shell.add_command('cd gw-xi')
    shell.add_command('rm -f pp.json pp/o-pp*') #cleanup
    shell.add_command('%s -F yambo_ppa.in -J pp -C pp' % yambo)
    shell.run()
    shell.clean()

    print ("Running Real Axis in folder 'gw-xi/ra'")
    ra = YamboIn('%s -d -g n -V all'%yambo,folder='gw-xi')
    ra['EXXRLvcs'] = [80,'Ry']       # Self-energy. Exchange
    ra['BndsRnXd'] = [1,25]          # Screening. Number of bands
    ra['NGsBlkXd'] = [3,'Ry']        # Cutoff Screening
    ra['GbndRnge'] = [1,25]          # Self-energy. Number of bands
    ra['QPkrange'][0][2:] = [2, 6]       # QP range. All BZ
    ra.write('gw-xi/yambo_ra.in')
    shell = bash() 
    shell.add_command('cd gw-xi')
    shell.add_command('rm -f ra.json ra/o-ra*') #cleanup
    shell.add_command('%s -F yambo_ra.in -J ra -C ra' % yambo)
    shell.run()
    shell.clean()

def plot_xi():
    #pack the files in .json files
    pack_files_in_folder('gw-xi')
    ya = YamboAnalyser('gw-xi')
    print('Plot Band structure for COHSEX, PPA and RA')
    path = [[[0,   0,   0],'$\Gamma$'],
            [[0.5, 0,   0],'M'],
            [[0.3333,0.3333, 0.0],'K'],
            [[0.0, 0.0, 0.0],'$\Gamma$']]
    ya.plot_gw_path('qp',path, cols=(lambda x: x[2]+x[3],))

def dyson_eq():
    #create the folder to run the calculation
    folder_dyson = 'gw-zeros'
    if not os.path.isdir(folder_dyson):
        shell = bash() 
        shell.add_command('mkdir -p %s' % folder_dyson)
        shell.add_command('cp -r database/SAVE %s/' % folder_dyson)
        shell.run()
        shell.clean()

    dyson = YamboIn('%s -p p -g n -V all'%yambo,folder=folder_dyson)

    dyson['EXXRLvcs'] = [80,'Ry']           # Self-energy. Exchange
    dyson['BndsRnXp'] = [1,25]              # Screening. Number of bands
    dyson['NGsBlkXp'] = [ 3,'Ry']        # Cutoff Screening
    dyson['GbndRnge'] = [1,25]              # Self-energy. Number of bands
    dyson['QPkrange'][0][2:] = [2, 6]

    dyson['DysSolver'] = "n" 
    dyson.write('%s/yambo_newton.in' % folder_dyson)
    dyson['DysSolver'] = "s" 
    dyson.write('%s/yambo_secant.in' % folder_dyson)
    shell = bash() 

    print("calculating with Newton and Secant solver in folder 'gw-zeros/newton' and 'gw-zeros/secant'...")
    shell.add_command('cd %s' % folder_dyson)
    shell.add_command('rm -f *.json newton/o-* secant/o-*') #cleanup
    shell.add_command('%s -F yambo_newton.in -J newton -C newton' % yambo)
    shell.add_command('%s -F yambo_secant.in -J secant -C secant' % yambo)
    shell.run()
    shell.clean()

def plot_dyson():
    #pack the files in .json files
    pack_files_in_folder('gw-zeros')
    ya = YamboAnalyser('gw-zeros')
    print('plot kpoints for Newton and secant solver')
    path = [[[0,   0,   0],'$\Gamma$'],
            [[0.5, 0,   0],'M'],
            [[0.3333,0.3333, 0.0],'K'],
            [[0.0, 0.0, 0.0],'$\Gamma$']]
    ya.plot_gw_path('qp',path, cols=(lambda x: x[2]+x[3],))

if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='GW convergence')
    parser.add_argument('-c'  ,'--convergence',  action="store_true", help='Run convergence calculations')
    parser.add_argument('-p'  ,'--plot', action="store_true",         help='Pack into json files and plot the convergence results')
    parser.add_argument('-g'  ,'--gw', action="store_true",           help='Run a single GW calculation')
    parser.add_argument('-r'  ,'--results', action="store_true",      help='Pack into json files and plot a single GW calculation')
    parser.add_argument('-x'  ,'--xi', action="store_true",           help='GW calculations for several approximations of the Screenning')
    parser.add_argument('-xp' ,'--xp', action="store_true",           help='Plot GW results for COHSEX, PPA and RA')
    parser.add_argument('-z'  ,'--zeros', action="store_true",        help='GW calculations for Newton and Secant Solver')
    parser.add_argument('-zp' ,'--zp', action="store_true",           help='Plot GW results for Newton and Secant Solver')

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    create_save()
    if args.convergence:    gw_convergence()
    if args.plot:           plot_convergence()
    if args.gw:             gw()
    if args.results:        plot_gw()
    if args.xi:             xi()
    if args.xp:             plot_xi()
    if args.zeros:          dyson_eq()
    if args.zp:             plot_dyson()
