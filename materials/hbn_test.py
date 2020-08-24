from qepy import *
from yambopy import *

class hBN_1l_test():
    """
    Input files for hBN monolayer for testing purposes (i.e., UNCONVERGED)

    """

    def __init__(self,prefix='bn_TEST',pseudo_path='./pseudos'):

        self.prefix = prefix
        self.pseudos = pseudo_path

        self.scfin = self.scf()
        self.nscfin = self.nscf()
        self.ipin = self.ip()

    def get_inputfile(self):
        """ Define a Quantum espresso input file for boron nitride
        """
        layer_separation = 12
        qe = PwIn()
        qe.set_atoms([['N',[0.0,0.0,0.5]],
                      ['B',[1/3,2/3,0.5]]])
        qe.atypes = {'B': [10.811, "B.pbe-mt_fhi.UPF"],
                     'N': [14.0067,"N.pbe-mt_fhi.UPF"]}

        qe.control['prefix'] = "'%s'"%self.prefix
        qe.control['verbosity'] = "'high'"
        qe.control['wf_collect'] = '.true.'
        qe.control['pseudo_dir'] = "'%s'"%self.pseudos
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

    #scf
    def scf(self):
        qe = self.get_inputfile()
        qe.control['calculation'] = "'scf'"
        qe.kpoints = [9,9,1]
        return qe

    #nscf
    def nscf(self):
        qe = self.get_inputfile()
        qe.control['calculation'] = "'nscf'"
        qe.electrons['diago_full_acc'] = ".true."
        qe.electrons['conv_thr'] = 1e-8
        qe.system['nbnd'] = 70
        qe.system['force_symmorphic'] = ".true."
        qe.kpoints = [18,18,1]
        return qe

    #ip
    def ip(self):
        yip = YamboIn()
        yip.arguments.append('chi')
        yip.arguments.append('dipoles') 
        yip.arguments.append('optics')
        yip['QpntsRXd']=[1,1]
        yip['DipBands']=[1,70]
        yip['NElectro']=8
        yip['DmRngeXd']=[[0.1,0.1],'eV']
        yip['ETStpsXd']=500
        yip['FFTGvecs']=[30,'Ry']
        yip['BndsRnXd']=[1,70]
        yip['EnRngeXd']=[[0.,10.],'eV']
        yip['LongDrXd']=[1.0,0.0,0.0]
        yip['Chimod']="IP" 
        return yip

    def __str__(self):
        s = ""
        s += "material: hBN monolayer for TESTING purposes\n"
        s += "input files available: qe-scf, qe-nscf, yambo-ip"
        return s
