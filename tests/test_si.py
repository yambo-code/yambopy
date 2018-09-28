from __future__ import print_function
#
# Author: Henrique Pereira Coutada Miranda
# Tests for yambopy
# Si
#
import matplotlib
import unittest
import sys
import os
import shutil
import argparse
import subprocess
import filecmp
import shutil as sh
from yambopy import *
from qepy import *

class TestPW_Si(unittest.TestCase):
    """ This class creates the input files for Si and compares them to reference files
    """
    def get_inputfile(self):
        qe = PwIn()
        qe.atoms = [['Si',[0.125,0.125,0.125]],
                    ['Si',[-.125,-.125,-.125]]]
        qe.atypes = {'Si': [28.086,"Si.pbe-mt_fhi.UPF"]}

        qe.control['prefix'] = "'si'"
        qe.control['wf_collect'] = '.true.'
        qe.control['pseudo_dir'] = "'../pseudos'"
        qe.system['celldm(1)'] = 10.3
        qe.system['ecutwfc'] = 40
        qe.system['occupations'] = "'fixed'"
        qe.system['nat'] = 2
        qe.system['ntyp'] = 1
        qe.system['ibrav'] = 2
        qe.kpoints = [4, 4, 4]
        qe.electrons['conv_thr'] = 1e-8
        return qe

    def test_pw_input_relax(self):
        """ Generate a silicon pw.x input file for the relaxation cycle
        """
        if not os.path.isdir('relax'):
            os.mkdir('relax')
        qe = self.get_inputfile()
        qe.control['calculation'] = "'vc-relax'"
        qe.ions['ion_dynamics']  = "'bfgs'"
        qe.cell['cell_dynamics']  = "'bfgs'"
        qe.write('relax/si.scf')
        self.assertEqual(filecmp.cmp('relax/si.scf', 'reference/si/relax_si.scf'),True)

    def test_pw_input_scf(self):
        """ Generate a silicon pw.x input file for the self consistent cycle
        """
        if not os.path.isdir('scf'):
            os.mkdir('scf')
        qe = self.get_inputfile()
        qe.control['calculation'] = "'scf'"
        qe.write('scf/si.scf')
        self.assertEqual(filecmp.cmp('scf/si.scf', 'reference/si/scf_si.scf'),True)

    def test_pw_input_nscf(self):
        """ Generate a silicon pw.x input file for the non self consistent cycle
        """
        if not os.path.isdir('nscf'):
            os.mkdir('nscf')
        qe = self.get_inputfile()
        qe.control['calculation'] = "'nscf'"
        qe.electrons['diago_full_acc'] = ".true."
        qe.electrons['conv_thr'] = 1e-8
        qe.system['nbnd'] = 30
        qe.system['force_symmorphic'] = ".true."
        qe.kpoints = [2, 2, 2]
        qe.write('nscf/si.nscf')
        self.assertEqual(filecmp.cmp('nscf/si.nscf', 'reference/si/nscf_si.nscf'),True)

class TestPW_Si_Run(unittest.TestCase):
    """ This class creates the input files and runs the pw.x code
    """
    def test_pw_si(sef):
        """ Run relaxation, self consistent cycle and non self consistent cycle
        """
        print("\nstep 1: relax")
        os.system('cd relax; pw.x < si.scf > si.scf.log')

        e = PwXML('si',path='relax')
        pos = e.get_scaled_positions()

        q = PwIn('scf/si.scf')
        print("old celldm(1)", q.system['celldm(1)'])
        q.system['celldm(1)'] = e.cell[0][2]*2
        print("new celldm(1)", q.system['celldm(1)'])
        q.atoms = list(zip([a[0] for a in q.atoms],pos))
        q.write('scf/si.scf')

        print("step 2: scf")
        os.system('cd scf; pw.x < si.scf > si.scf.log')
        os.system('cp -r scf/si.save nscf')

        print("step 3: nscf")
        os.system('cd nscf; pw.x < si.nscf > si.nscf.log')


class TestYamboPrep_Si(unittest.TestCase):
    def test_yambo_preparation(self):
        """ Run p2y and yambo to prepare the database
        """
        if not os.path.isdir('database'):
            os.mkdir('database')
            os.system('cd nscf/si.save; p2y 2> ../../database/p2y.log')
            os.system('cd nscf/si.save; yambo 2> ../../database/yambo.log')
            os.system('mv nscf/si.save/SAVE database')

class TestYamboIn_GW_Si(unittest.TestCase):
    def setUp(self):
        """ Prepare the databases
        """
        if not os.path.isdir('database/SAVE'):
            os.makedirs('database')
            os.system('cd database; tar xfz ../reference/si/yambo_gw_conv/gw_conv.tar.gz')
        if not os.path.isdir('gw_conv/SAVE'):
            sh.copytree('database/SAVE','gw_conv/SAVE')
        if not os.path.isdir('gw/SAVE'):
            sh.copytree('database/SAVE','gw/SAVE')

    def test_gw_input(self):
        """ Test if we can initialize the YamboIn class for a typical GW input file
        """
        y = YamboIn('yambo -p p -g n -V all',folder='gw')

    def test_gw_convergence(self):
        """ Test if we can generate multiple input files changing some variables
        """
        y = YamboIn('yambo -p p -g n -V all',folder='gw_conv')
        conv = { 'FFTGvecs': [[5,10,15],'Ry'],
                 'NGsBlkXp': [[1,2,5], 'Ry'],
                 'BndsRnXp': [[1,10],[1,20],[1,30]] }
        y.optimize(conv)
        return y


class TestYamboIn_GW_Si_Run(unittest.TestCase):
    def test_yambo_gw_si(self):
        """ Run GW calculation with yambo
        """
        y = YamboIn('yambo -p p -g n -V all',folder='gw_conv')
        conv = { 'FFTGvecs': [[5,10,15],'Ry'],
                 'NGsBlkXp': [[1,2,5], 'Ry'],
                 'BndsRnXp': [[1,10],[1,20],[1,30]] }
        y.optimize(conv)

        def run(filename):
            folder = filename.split('.')[0]
            print(filename, folder)
            os.system('cd gw_conv; yambo -F %s -J %s -C %s 2> %s.log'%(filename,folder,folder,folder))

        y.optimize(conv,run=run)

    def test_yambopy_analysegw(self):
        """ Test the yambopy analysegw executable
        """
        os.system('yambopy analysegw gw_conv FFTGvecs -bc 5 -kc 3 -bv 4 -kv 1 -nd')
        out = np.loadtxt('analyse_gw_conv/gw_conv_FFTGvecs.dat')
        ref = np.loadtxt('reference/si/analyse_gw_conv/gw_conv_FFTGvecs.dat')
        print("ref:")
        print(ref)
        print("out:")
        print(out)
        self.assertEqual(np.isclose(ref,out,atol=1e-3).all(),True)

        os.system('yambopy analysegw gw_conv BndsRnXp -bc 5 -kc 3 -bv 4 -kv 1 -nd')
        out = np.loadtxt('analyse_gw_conv/gw_conv_BndsRnXp.dat')
        ref = np.loadtxt('reference/si/analyse_gw_conv/gw_conv_BndsRnXp.dat')
        print("ref:")
        print(ref)
        print("out:")
        print(out)
        self.assertEqual(np.isclose(ref,out,atol=1e-3).all(),True)


class TestYamboIn_BSE_Si(unittest.TestCase):
    def setUp(self):
        """ Prepare the databases
        """
        if not os.path.isdir('database/SAVE'):
            os.makedirs('database')
            os.system('cd database; tar xfz ../reference/si/yambo_bse_conv/bse_conv.tar.gz')
        if not os.path.isdir('bse/SAVE'):
            sh.copytree('database/SAVE','bse/SAVE')
        if not os.path.isdir('bse_conv/SAVE'):
            sh.copytree('database/SAVE','bse_conv/SAVE')

    def test_bse_input(self):
        """ Test if we can initialize the YamboIn class for a typical BSE input file
        """
        y = YamboIn('yambo -b -o b -k sex -y h -V all',folder='bse')

    def test_bse_convergence(self):
        """ Test if we can generate multiple input files changing some variables
        """
        y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse_conv')
        conv = { 'FFTGvecs': [[5,10,15],'Ry'],
                 'NGsBlkXs': [[1,2,5], 'Ry'],
                 'BndsRnXs': [[1,10],[1,20],[1,30]] }
        y.optimize(conv)
        return y

class TestYamboIn_BSE_Si_Run(unittest.TestCase):
    def test_yambo_bse_si(self):
        """ Run BSE calculation with yambo
        """
        y = YamboIn('yambo -b -o b -k sex -y d -V all',folder='bse_conv')
        conv = { 'FFTGvecs': [[5,10,15],'Ry'],
                 'NGsBlkXs': [[1,2,5], 'Ry'],
                 'BndsRnXs': [[1,10],[1,20],[1,30]] }

        print()
        def run(filename):
            folder = filename.split('.')[0]
            print(filename, folder)
            os.system('cd bse_conv; yambo -F %s -J %s -C %s 2> %s.log'%(filename,folder,folder,folder))

        y.optimize(conv,run=run)

class TestYamboOut_BSE_Si(unittest.TestCase):
    def test_yamboout_bse_si(self):
        """ Read the yambo BSE output files and write them as .json
        """
        for dirpath,dirnames,filenames in os.walk('bse_conv'):
            #check if there are some output files in the folder
            if ([ f for f in filenames if 'o-' in f ]):
                y = YamboOut(dirpath,save_folder='bse_conv')
                y.pack()

    def test_yamboanalyse_bse_si(self):
        """ Analyse the BSE .json output files
        """
        y = YamboAnalyser('bse_conv')
        y.plot_bse('eps')

    def test_yambopy_analysebse(self):
        """ Test the yambopy analysebse executable
        """
        os.system('yambopy analysebse bse_conv FFTGvecs -nd')
        out = np.loadtxt('analyse_bse_conv/bse_conv_FFTGvecs_excitons.dat')
        ref = np.loadtxt('reference/si/analyse_bse_conv/bse_conv_FFTGvecs_excitons.dat')
        print("ref:")
        print(ref)
        print("out:")
        print(out)
        self.assertEqual(np.isclose(ref,out,atol=1e-3).all(),True)

        os.system('yambopy analysebse bse_conv BndsRnXs -nd')
        out = np.loadtxt('analyse_bse_conv/bse_conv_BndsRnXs_excitons.dat')
        ref = np.loadtxt('reference/si/analyse_bse_conv/bse_conv_BndsRnXs_excitons.dat') 
        print("ref:")
        print(ref)
        print("out:")
        print(out)
        self.assertEqual(np.isclose(ref,out,atol=1e-3).all(),True)

class TestYamboOut_GW_Si(unittest.TestCase):
    def test_yamboout_gw_si(self):
        """ Read the yambo GW output files
        """
        for dirpath,dirnames,filenames in os.walk('gw_conv'):
            #check if there are some output files in the folder
            if ([ f for f in filenames if 'o-' in f ]):
                y = YamboOut(dirpath,save_folder='gw_conv')
                y.pack()

    def test_yamboanalyse_gw_si(self):
        """ Analyse the yambo GW .json output files
        """
        y = YamboAnalyser('gw_conv')
        y.plot_gw('qp')

if __name__ == '__main__':
    #parse options
    parser = argparse.ArgumentParser(description='Test the yambopy script.')
    parser.add_argument('-i','--input', action="store_true",
                        help='Generate the input files and compare with the reference ones')
    parser.add_argument('-f','--full',  action="store_true",
                        help='Generate the input files, run them and compare the results')
    parser.add_argument('-c','--clean',  action="store_true",
                        help='Clean all the data from a previous run')
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    #first test if yambo is installed
    sp = subprocess.PIPE
    yambo_not_available = subprocess.call("yambo", shell=True, stdout=sp, stderr=sp)
    if yambo_not_available:
        print("yambo not found, please install it before running the tests")
        sys.exit(1)

    # Count the number of errors
    nerrors = 0
    ul = unittest.TestLoader()
    tr = unittest.TextTestRunner(verbosity=2)

    #
    # Test pw.x
    #
    suite = ul.loadTestsFromTestCase(TestPW_Si)
    nerrors += not tr.run(suite).wasSuccessful()

    if args.full:
        suite = ul.loadTestsFromTestCase(TestPW_Si_Run)
        nerrors += not tr.run(suite).wasSuccessful()

    #
    # Test p2y and yambo
    #
    if args.full:
        suite = ul.loadTestsFromTestCase(TestYamboPrep_Si)
        nerrors += not tr.run(suite).wasSuccessful()

    #
    # Test GW on yambo
    #
    suite = ul.loadTestsFromTestCase(TestYamboIn_GW_Si)
    nerrors += not tr.run(suite).wasSuccessful()

    if args.full:
        suite = ul.loadTestsFromTestCase(TestYamboIn_GW_Si_Run)
        nerrors += not tr.run(suite).wasSuccessful()

        suite = ul.loadTestsFromTestCase(TestYamboOut_GW_Si)
        nerrors += not tr.run(suite).wasSuccessful()

    #
    # Test BSE on yambo
    #
    suite = ul.loadTestsFromTestCase(TestYamboIn_BSE_Si)
    nerrors += not tr.run(suite).wasSuccessful()

    if args.full:
        suite = ul.loadTestsFromTestCase(TestYamboIn_BSE_Si_Run)
        nerrors += not tr.run(suite).wasSuccessful()

        suite = ul.loadTestsFromTestCase(TestYamboOut_BSE_Si)
        nerrors += not tr.run(suite).wasSuccessful()

    #clean tests
    if args.clean or nerrors==0:
        print("cleaning...")
        os.system('rm -rf scf bse bse_conv gw gw_conv nscf relax database '
                  'analyse_bse_conv analyse_gw_conv proj.in')
        print("done!")

    sys.exit(nerrors)
