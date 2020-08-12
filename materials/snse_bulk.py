from qepy import *
from yambopy import *

class SnSe_bulk():
    """
    Input files for SnSe bulk
    
    Notes:
        - Pseudopotentials from pseudo dojo
                (PBE, norm-conserving, full relativistic, standard accuracy)
        - Atomic positions at EXPERIMENTAL coordinates [Sint & al., Act. Cryst. B (2016)]
        - Pcmn lattice axes orientation is used
        - Spin-orbit is omitted [but can be easily added]
    """

    def __init__(self,prefix='snse_nosoc',pseudo_path='./pseudos'):

        self.prefix = prefix
        self.pseudos = pseudo_path

        self.scf = scf()
        self.nscf = nscf()
        self.ip = ip()

    def get_inputfile(self):
        """ Define a Quantum espresso input file for tin selenide
        """
        qe = PwIn()
        qe.set_atoms([['Sn',[0.89703,0.75,0.88151]],
                      ['Sn',[0.39703,0.25,0.61849]],
                      ['Sn',[0.10297,0.25,0.11849]],
                      ['Sn',[0.60297,0.75,0.38151]],
                      ['Se',[0.51796,0.75,0.14477]],
                      ['Se',[0.01796,0.25,0.35523]],
                      ['Se',[0.48204,0.25,0.85523]],
                      ['Se',[0.98204,0.75,0.64477]]])
        qe.atypes = {'Sn': [118.71, "Sn.pbe_nc_fr_standard_dojo.upf"],
                     'Se': [78.96,"Se.pbe_nc_fr_standard_dojo.upf"]}

        qe.control['prefix'] = "'%s'"%self.prefix
        qe.control['verbosity'] = "'high'"
        qe.control['wf_collect'] = '.true.'
        qe.control['pseudo_dir'] = "'%s'"%self.pseudos
        qe.system['celldm(1)'] = 8.3936904098
        qe.system['celldm(2)'] = 0.9345325604
        qe.system['celldm(3)'] = 2.5877570777
        qe.system['ecutwfc'] = 80
        qe.system['occupations'] = "'fixed'"
        qe.system['vdw_corr']= "'grimme-d2'"
        qe.system['nat'] = 8
        qe.system['ntyp'] = 2
        qe.system['ibrav'] = 8
        qe.electrons['conv_thr'] = 1e-8

        qe.kpoints = [9,9,3]
        return qe

    #scf
    def scf(self):
        qe = self.get_inputfile()
        qe.control['calculation'] = "'scf'"
        qe.kpoints = [9,9,3]
        return qe

    #nscf
    def nscf(self):
        qe = self.get_inputfile()
        qe.control['calculation'] = "'nscf'"
        qe.electrons['diago_full_acc'] = ".true."
        qe.electrons['conv_thr'] = 1e-8
        qe.system['nbnd'] = 70
        qe.system['force_symmorphic'] = ".true."
        qe.kpoints = [18,18,6]
        return qe

    #ip
    def ip(self):
        yip = YamboIn()
        yip.arguments.append('chi')
        yip.arguments.append('dipoles') 
        yip.arguments.append('optics')
        yip['DIP_ROLEs']= "'k c v'"
        yip['DIP_CPU'] = "'32 1 2'"
        yip['X_ROLEs']= "'q g k c v'"
        yip['X_CPU'] = "'1 1 32 1 2'"
        yip['QpntsRXd']=[1,1]
        yip['DipBands']=[50,70]
        yip['DmRngeXd']=[[0.05,0.05],'eV']
        yip['ETStpsXd']=2000
        yip['FFTGvecs']=[30,'Ry']
        yip['BndsRnXd']=[50,70]
        yip['EnRngeXd']=[[0.,3.],'eV']
        yip['LongDrXd']=[0.61237,0.61237,0.5] #Pump with laser direction at 30 degs with normal to the sample 
    
        #     \  |      #####################
        #      \ |      #  sqrt(6)/4 ^  
        #       \|      #            |__ >
        ##############  #
        ### SAMPLE ###  # SAMPLE PLANE xy
        ##############  #####################

        yip['Chimod']="'IP'" 
        yip['XfnQP_E']=[0.368058,1.,1.] #This matches the fundamental GW gap with that of Shi and Kioupakis
        #yip['XfnQP_E']=[0.47732,1.,1.] #This matches the minimum direct GW gap with that of Shi and Kioupakis
        return yip

    def __str__(self):
        s = ""
        s += "material: SnSe bulk\n"
        s += "input files available: qe-scf, qe-nscf, yambo-ip"
        return s
