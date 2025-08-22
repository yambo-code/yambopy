#import libs 
import os 
import sys
from netCDF4 import Dataset
import numpy as np 
from collections import defaultdict
from tqdm import tqdm 
import h5py
from yambopy.pert2yambo.grid_match import GridMatch 
from itertools import product 
import pickle 
import yaml  

#read perturbo data and account for missing occupations

class dynamic_occupations():
  
    #store indices and number of bands ndb
  
    vb_idx = None
    cb_idex = None
    occup_template = None
    
    full_num_bands = None
    
    #get info on how mnay pert kpts there are in the scattering channels
    
    kpts_used = None
    
    #store the band max and min 
    
    pert_band_min = None 
    pert_band_max = None 
    
    def __init__(self,teth5file,cdynafile,ndbfile,dyn_yamlfile,tmp_out="./tmp"):
        
        self.teth5file = teth5file
        self.cdynafile = cdynafile
        self.tmp_out = tmp_out
        self.ndbfile = ndbfile
        self.dyn_yamlfile = dyn_yamlfile
    
    #def get_temper_data(self):
        
        #f = open(self.temperfile)
        
        #chars = f.readlines()[1].split()
        #self.mu = float(chars[1])
        #self.temp = float(chars[0])
    
    #@staticmethod
    #def fermi_dirac(E, mu=0.0, T=300.0, kB=8.617333262145e-5):
    #    beta = 1 / (kB * T)
    #    return 1.0 / (1.0 + np.exp(beta * (E - mu)))
    
    
    #def read_zero_temp_data(self):
        #"""
        
        #Retrieves zero_temp DFT energies and occupations from the yambo ndbcarriers file. 
        
        #"""
        
        #with Dataset(self.ndbfile) as ndbfile:
        #    zero_temp_energies = ndbfile.variables["RT_carriers_E_bare"][()]
        #    zero_temp_occupations = ndbfile.variables["RT_carriers_f_bare"][()]
        #    
        #    num_kpts,_,num_bands = ndbfile.variables['RT_carriers_dimensions'][()][0:3]

            #reshape the data so we get all zero temps energies as kpt: E_0, . . . E_N
        #    zero_temp_energies.reshape(num_kpts,num_bands)
        #    zero_temp_occupations.reshape(num_kpts,num_bands)
            
        #    return zero_temp_energies,zero_temp_occupations,num_kpts,num_bands

    #def read_ndb_kpts(self):
    #    """
    #    Get full grid of points from the ndb file, and returns the ordered array of kpts. 
        
    #    """
    #    with Dataset(self.ndbfile) as ndbfile:
    #        return ndbfile.variables['RT_kpt'][()].T
    
            
    #def create_ndb_dict(self):
        
    #    """
        
    #    Creates a dict of dicts with structure: 
        
    #    kpt:{E:[],f:[]}
        
    #    """
    #    self.ndb_data = defaultdict(dict)
        
    #    zero_temp_energies,occupations,num_kpts,num_bands = self.read_zero_temp_data()
        
    #    #calculate equilibrium occupations with the Fermi-Dirac distribution
    #    occupations = self.fermi_dirac(zero_temp_energies,self.mu,self.temp)
        
    #    #create dictionary of dictrionaries with np.array entries for the respective quantities
    #    for kpt,E,f in zip(self.read_ndb_kpts(),zero_temp_energies,occupations):
            
    #        key = tuple(kpt.tolist())
            
    #        self.ndb_data[key]["E"] = E
    #        self.ndb_data[key]["f"] = f
        
        
    def read_perturbo_kpts(self):
        """
        
        Read kpts from perturbo prefix_tet.kpt file. If the kpts were reduced by perturbo, then it will creake two .pkl files
        accoutning for the difference; else, it simply reads the entire BZ from perturbo. 
        
        """
        with h5py.File(self.teth5file, 'r') as teth5:
            
            return teth5["kpts_all_crys_coord"][:,:]    
            
    def read_pert_num_kpts(self):
        
        """
        Get number of num of perturbo kpts. 
       
        """
        
        with h5py.File(self.teth5file, 'r') as teth5:
            self.num_pert_kpts = teth5['num_kpts'][()]
    
    
    #def update_dynamic_entries(self):
    
    def get_vcb_indices(self):
        """
        Returns the indices of the valence and conduction bands.
         
        """
            
        with Dataset(self.ndbfile) as ndbfile:
            
            num_kpts,_,num_bands = ndbfile.variables['RT_carriers_dimensions'][()][0:3]
            occupations = ndbfile.variables['RT_carriers_f_bare'][()].reshape(num_kpts,num_bands)

            #assign valence band index, and conduction band index
            self.vb_idx = max(np.where(occupations[0] >= 1.0)[0])
            self.cb_idx = min(np.where(occupations[0] == 0.0)[0])
            
            self.full_num_bands = num_bands   
            self.num_kpts = num_kpts

            print("Number of Yambo k-points ",num_kpts)
            print("Number of Yambo bands    ",num_bands)
            print("Top valence:    ",self.vb_idx+1)
            print("Bottom conduction:     ",self.cb_idx+1)
            
            #construct occupation template for later, and assign occups as one would normally see in QE
            self.occup_template = occupations[0]
                
                
            


    def parse_bands_from_yaml(self):
        with open(self.dyn_yamlfile, 'r') as f:
            data = yaml.safe_load(f)
            
            input_params = data.get('input parameters', {})
            before_conv = input_params.get('before conversion (before reading from epr file)', {})

            self.pert_band_min = before_conv.get('band_min')
            self.pert_band_max = before_conv.get('band_max')

            self.pert_bands = list(range(self.pert_band_min,self.pert_band_max+1))
    
    def pert_grid_reduced(self):
        
        with h5py.File(self.teth5file, 'r') as teth5:
            self.kpts_used = teth5['num_kpts'][()]
            self.kpts_grid = teth5['kgrid_dim'][()]
            self.kpts_grid_size = np.prod(self.kpts_grid)
                        
            return True if self.kpts_used != self.kpts_grid_size else False
    
    @staticmethod
    def create_grid(kdims):
        """
        Generator-based Monkhorst-Pack grid for perturbo calculations.
        Yields each k-point one at a time instead of storing in memory.
        """
        grids = (np.linspace(0.0, 1.0, kdim, endpoint=False) for kdim in kdims)
        return product(*grids)  
    
    
    def get_files(self,time_range=[],run=1,shift=[0,0,0]):
            """
            
            """
            
            #create tmp folder
            if not os.path.exists(self.tmp_out):
                os.mkdir(self.tmp_out)
            
            #check if the perturbo grid is reduced
            
            reduced = self.pert_grid_reduced()
            
            with h5py.File(self.cdynafile, "r") as cdyna:
                
                #check time range if empty
                if len(time_range) == 0:
                    
                    time_range.append(0)
                    time_range.append(len(cdyna[f'dynamics_run_{run}'].keys())-2)

                #loop over all snapshots
                self.occupation = {}
                print("Reading occupations...")
                for snapshot_t in tqdm(cdyna[f'dynamics_run_{run}'].keys()):
                    
                    if snapshot_t == 'num_steps' or snapshot_t == 'time_step_fs':
                        continue
                    
                    #get filenames to store the occupations
                    occupsfilename = os.path.join(self.tmp_out,f"dynamic_occups_run_{1}_snap_{snapshot_t}.pkl")
#                    print("occups file name ",occupsfilename)
                    
                    if reduced is True:
                        global red_occupsfilename
                        red_occupsfilename = os.path.join(self.tmp_out,f"dynamic_occups_run_{1}_snap_{snapshot_t}_redkpts.pkl")
                    
                    indx=int(snapshot_t.split('_')[-1])
#                    print("Show occupations :",str(snapshot_t))
#                    print(cdyna[f'dynamics_run_{run}'][snapshot_t][:,:])
                    #with open(red_occupsfilename, "ab") as red_occupspkl:
                    self.occupation[indx]=cdyna[f'dynamics_run_{run}'][snapshot_t][:,:]
                    
                    with open(occupsfilename,"wb") as occupspkl:
                        snap = int(snapshot_t.split("_")[-1])
                        
                        #verify that the occupations are not zero
                        pert_sum_occups = np.sum(cdyna[f'dynamics_run_{run}'][snapshot_t][:,:])
#                        assert pert_sum_occups > 0.0, "Bands have occupation of zero!" #need to verify that the bands are somewhat occupied
                        
                                            
                        #loop over relevant timesteps
                        if snap >= min(time_range) and snap <= max(time_range):

                            #define occupations dict for pert points --> used later
                        
                            pert_occups = dict()                        

                            #open .pkl for the reduced grid
                            if reduced is True:
                                global red_occupspkl
                                red_occupspkl = open(red_occupsfilename,"wb")

                            #retrieve the perturbo occupations and create a dict; write if reduced
#                            with h5py.File(self.teth5file, 'r') as teth5:
#                                for index in range(0,self.kpts_used):
                                
#                                    kpt = teth5["kpts_all_crys_coord"][index,:]
#                                    key = tuple(kpt.tolist())
#                                    
#                                    occups = cdyna[f'dynamics_run_{run}'][snapshot_t][index,:]
#                                    pert_occups[key] = occups
                                
                            
#                                    if reduced is True:
#                                        pickle.dump([kpt,occups],red_occupspkl)
                            
#                            #close file
#                            if reduced is True:
#                                red_occupspkl.close()  
                                       
                            #assign the occupations to the full perturbo grid and write to file
#                            with open(occupsfilename,"wb") as occupspkl:
#                            
#                                for full_kpt in tqdm(self.create_grid(self.kpts_grid)):
#                                    
#                                    occups = self.occup_template #assign baseline occupations at zero temp
#                                    
#                                    try:
#                                        occups = pert_occups[full_kpt]
#                                            
#                                    except KeyError:
#                                        pass
#                                    
#                                    pickle.dump([full_kpt,occups],occupspkl)
                
                                    
    

if __name__ == "__main__":
    
    ndb = "/Users/dani/projects/GaAS_workflow/calcs/ver_calcs/080525_occupation_read_test/ndb.RT_carriers"

    cdyna = "/Users/dani/projects/GaAS_workflow/calcs/ver_calcs/200625_yambopy_RTBD_test/gaas_cdyna.h5"
    teth5 = "/Users/dani/projects/GaAS_workflow/calcs/ver_calcs/200625_yambopy_RTBD_test/gaas_tet.h5"
    
    yamlfile = "/Users/dani/projects/GaAS_workflow/calcs/ver_calcs/200625_yambopy_RTBD_test/gaas_dynamics-run.yml" 
    
    dynoccups = dynamic_occupations(tmp_out="./tmp",dyn_yamlfile=yamlfile,cdynafile=cdyna,teth5file=teth5,ndbfile=ndb)
    
    dynoccups.get_vcb_indices()
    dynoccups.parse_bands_from_yaml()
    dynoccups.pert_grid_reduced()
    #print(dynoccups.kpts_grid)
    dynoccups.get_files()
