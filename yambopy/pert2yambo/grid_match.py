# import libs 
import numpy as np
from collections import defaultdict, Counter
from tqdm import tqdm 
import csv, os, pickle
import matplotlib.pyplot as plt
#create class for reversible dictionary

class ReversibleDict(dict):
    def __setitem__(self, key, value):
        # Remove any previous connections with these values
        if key in self:
            del self[key]
        if value in self:
            del self[value]
        dict.__setitem__(self, key, value)
        dict.__setitem__(self, value, key)

    def __delitem__(self, key):
        dict.__delitem__(self, self[key])
        dict.__delitem__(self, key)

    def __len__(self):
        """Returns the number of connections"""
        return dict.__len__(self) // 2
    
    def max_by_position(self):
        """Returns a tuple of maximum values per position across all unique tuples (keys and values)."""
        seen = set()
        elements = []

        for k, v in self.items():
            if k in seen or v in seen:
                continue
            seen.add(k)
            seen.add(v)
            elements.extend([k, v])

        if not elements:
            raise ValueError("ReversibleDict is empty.")

        # Transpose and compute max per position
        transposed = zip(*elements)
        return tuple(max(group) for group in transposed)


#define custom product to return numpy arrays 

def product(*iterables, repeat=1):
    """
    Customized function which returns all permutated Cartesian products of a a set of arrays.
    
    """
    if repeat < 0:
        raise ValueError('repeat argument cannot be negative')
    pools = [tuple(pool) for pool in iterables] * repeat

    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]

    for prod in result:
        yield np.array(prod)

#custom displacement generator function

def generate_displacements(points_between):
    """
    Given the number of dense points between coarse grid points per dimension,
    return symmetric displacements centered at 0, *including boundary points*.

    Args:
        points_between (tuple of ints): Number of dense points between coarse
        grid points along each dimension.

    Returns:
        np.ndarray: Array of displacements centered at zero.
    """
    ranges = []
    for n in points_between:
        # Add 1 to ensure we include boundary point symmetrically
        span = n if n % 2 == 0 else n + 1
        half_span = span // 2
        ranges.append(range(-half_span, half_span + 1))  # includes both ends
    return np.array(list(product(*ranges)))

#custom wrapping function for periodic indices later on 

def wrap_to_range(c, a, b):
    """
    c --> number that needs mapping
    a --> start of periodic range
    b --> end of periodic range
    
    """
    width = b - a
    return a + ((c - a) % width)

#create class for grid
class GridMatch():

    #match_dict = None
    dense_usage_count = Counter()
    num_lines = 0
    tmp_neighbours = None
    
    min_dist_array_1 = None
    kpt_different_1 = None
    min_dist_array_2 = None
    kpt_different_2 = None

    def find_lattice_info(self,lattice):
        """
        
        Find the smallest distance between kpoints and all the different k-points along each direction. 
        
        """
        
        # Find unique values along each axis
        kpt_different = defaultdict(set)
        for point in lattice:
            for i in range(3):
                kpt_different[i].add(point[i])

        kpt_different_array = np.array([len(kpt_different[i]) for i in range(3)])
        
        min_dist_array = 1.0/kpt_different_array
        return min_dist_array, np.array(kpt_different_array)
        
        
            
    def create_int_map(self,lattice,min_dist_array):
        """
        
        Map each point of the lattice to the nearest integer by dividing by the smallest distance, using a reversible dictionary. 
        
        kpt <--> integer multiple
        
        """      
        
        int_map = ReversibleDict()
         
        for kpt in tqdm(lattice):
            int_map[tuple(kpt.tolist())] = tuple(np.rint(kpt/min_dist_array).tolist())
        
        return int_map
        
    
    def get_nearest_neighbours(self,lattice1,lattice2,tmp="./tmp",filename="neighbours.csv",symm=True):
        """
        
        Given a coarse lattice 1, and a dense lattice 2, find the nearest neighbour in lattice 1 to each point in lattice 2, as well as the weight of each dense point with respect its coarse neighbours.
        
        Assumes that lattice1.shape > lattice2.shape
        
        Returns:
        
        Has a memory option to store all lines in a .csv file, useful for large grids. 
        
        Symm option allows to track the wrapping symmetry of the grid. 
        
        """
        
        assert lattice1.shape[0] < lattice2.shape[0] #check if lattice 1 is coarser than lattice 2 
        
        print(40*"-")
        
        #create int maps
        print("\n Retrieving information for coarse lattice . . . \n")
        self.min_dist_array_1,self.kpt_different_1 = self.find_lattice_info(lattice1) #coarse lattice
        int_map_1 = self.create_int_map(lattice1 ,self.min_dist_array_1)
        print(f"\n## Minimum distance between points: {self.min_dist_array_1} ## \n## Kpts in each direction: {self.kpt_different_1} ## \n## Kpts in total: {len(lattice1)} ##")
        
        print(40*"-")
        
        print("\n Retrieving information for dense lattice . . .\n")
        self.min_dist_array_2,self.kpt_different_2 = self.find_lattice_info(lattice2) #dense lattice
        int_map_2 = self.create_int_map(lattice2 ,self.min_dist_array_2)  #necessary since this contains also the coarse lattice points
        print(f"\n## Minimum distance between points: {self.min_dist_array_2} ## \n## Kpts in each direction: {self.kpt_different_2} ## \n## Kpts in total: {len(lattice2)} ##")
        
        print(40*"-")
        
        #find number of neighbours in each direction
        points_from_coarse = (self.kpt_different_2 // self.kpt_different_1).astype(int) + 1 #np.array([i/j for i,j in zip(kpt_different_2,kpt_different_1)])   
        
        #generate all displacements
        disps = generate_displacements(points_from_coarse)
        
        print(f"\n Neighbours in each direction: {points_from_coarse} \n")
        
        #find nearest neighbours by density; map the dense point to the coarse point, and include the weight
        
        print(40*"-")
        
        #create tmp folder
        if not os.path.exists(tmp):
            os.mkdir(tmp)
        
        print(f"\n Created tmp folder at {tmp} \n")
        
        # write the headers and create tmp file to store neighbours
        self.tmp_neighbours = os.path.join(tmp,"neighbours.pkl")
        
        if os.path.exists(self.tmp_neighbours):
            print(60*"#" + f"\n\nWarning: previous .pkl file found at {self.tmp_neighbours} . Removing entry . . .\n\n" + 60*"#")
            os.remove(self.tmp_neighbours)
        
        print(40*"-")
        
        print(f"\n Finding nearest neighbours . . . \n")
        

        #find nearest neighbours
        self.num_lines = 0
        
        with open(self.tmp_neighbours,"ab") as neighbourfile:
            for coarsept in tqdm(lattice1):
                
                #get dense integer representation of the coarsepoint 
                coarse_int_array = np.array(int_map_2[tuple(coarsept.tolist())])
                
                for displacement in disps:

                    neighbour = coarse_int_array + displacement #get the integer distance from the coarse center
        
                    #catch negative kpts and convert them to corresponding periodic representation
                    if symm is True:
                        if np.any(neighbour < 0) or np.any(neighbour >= self.kpt_different_2):
                            for index in np.where(neighbour < 0)[0]:
                                neighbour[index] = wrap_to_range(neighbour[index],0,self.kpt_different_2[index])#max_pos[index]) #maps the index out of range to one in range, as k-space is periodic
                            
                            for index in np.where(neighbour >= self.kpt_different_2)[0]:
                                neighbour[index] = wrap_to_range(neighbour[index],0,self.kpt_different_2[index])#,max_pos[index]) #maps the index out of range to one in range, as k-space is periodic
                            
                    #check if dense kpt neighbour exists --> otherwise not allowed 
                    try:
                        global dense_kpt
                        dense_kpt = int_map_2[tuple(neighbour.tolist())]
                    
                    
                    except KeyError:
                        continue
                    
                    #write to .pkl file
                    
                    pickle.dump([coarsept,dense_kpt],neighbourfile)
                    self.num_lines += 1 #useful for later
                    
                    #update histogram of points
                    self.dense_usage_count[dense_kpt] += 1
                    
                    
                    
            
        print("\n" + 40*"-")
        
        print(f"Assigning weights . . .")
        
        print(40*"-")
    
        #import pdb; pdb.set_trace()
    
    
        #keep track of weights
        weight_sum = 0
        count = 0 
        
        with open(self.tmp_neighbours, 'rb') as infile, open(filename, 'w') as outfile:
            writer = csv.writer(outfile)
            
            #write final header
            writer.writerow(["k1","k2","k3","k1'","k2'","k3'","weight"])
                
            for _ in tqdm(range(0,self.num_lines)):
                try:
                    data = pickle.load(infile)
                    dense_kpt = data[1]
                    coarse_kpt = data[0]
                
                    assert self.dense_usage_count[dense_kpt] > 0, f"Point {dense_kpt} does not appear!"
                
                    weight = 1.0/self.dense_usage_count[dense_kpt]
                    weight_sum += weight #update for weight tracking
                    count += 1
                    
                    assert weight <= 1.0, f"Weight cannot be more than unity! Got {weight}"
                    
                    #write the weights and neighbours to file
                    writer.writerow([*coarse_kpt,*dense_kpt,weight])
                    
                    
                except EOFError:
                    break
                
        assert count == self.num_lines, f"Not all lines of the .pkl were read! Red only {count} lines . . ."

        #import pdb;pdb.set_trace()

        print(40*"-")
        
        print(f"\n Nearest neighbours written to {filename} . . . \n")

        assert abs(weight_sum - len(lattice2)) < 1e-6, f"Normalization mismatch! Got {weight_sum} vs {len(lattice2)}"
        
        print(f"\nNormalization correct!\n")
        
        print("\n" + 70*"#")
        print("\n" + 30*"#" + " DONE " + 30*"#")
        print("\n" + 70*"#")
    
    def gen_coarse_occupations(self,dense_occupation_dict,band_labels=None, neighbourfilename=None, occupationfilename="coarse_occupations.csv"):
        """
        
        Generates occupations for the coarse points, given the occupations of the dense points. Assumes normalization of occupations
        in dense_occupation_dict.
        
        dense_occupation dict needs to be a dict of dicts of the form: {kpt_n:E_0, E_1 . . . E_m} 
        
        Where kpt_n is the tuple. 
        
        """

        with open(occupationfilename,"w") as occupationfile:
            
            
            print(40*"-")
        
            print(f"\nReading neighbours from {neighbourfilename} . . .\n")
        
            print(40*"-")

            #get the number of bands used
            
            num_bands = len(next(iter(dense_occupation_dict.values()))) 

            #extract info from neighbour files

            occups_dict = defaultdict(lambda: [0.0] * num_bands)

            print(40*"-")
        
            print("\nGenerating occupations for coarse points . . .\n")
        
            print(40*"-")

            #determine occupations and hold them in memory
            total_weighted_occupation = 0.0 #used to check normalization later on
            
            #import pdb; pdb.set_trace()
            
            if self.tmp_neighbours == None:
                neighbourfile = open(neighbourfilename,"r")
                
                reader = csv.reader(neighbourfile,delimiter=",")
                next(reader) # skip header
                
                for ln in reader:
                
                    ln = np.genfromtxt(ln)
                    
                    coarsept,densept,weight = ln[0:3],ln[3:6],ln[6:7]
                    
                    occup = weight[0] * dense_occupation_dict[tuple(densept)]
                    
                    occups_dict[tuple(coarsept)] += occup
                    total_weighted_occupation += occup #to keep track of weighted occupation numbers

                #NEED TO MIRROR THE FUNCTIONALITY BELOW 

            
            #load from pickle file for faster performance
            
            else:
                print("Found .pkl file from nearest neighbours . . .")
                
                assert self.tmp_neighbours.endswith(".pkl")
                with open(self.tmp_neighbours, 'rb') as pkfile:
                    for _ in tqdm(range(self.num_lines)):
                        try:
                            data = pickle.load(pkfile)
                            densept = data[1]
                            coarsept = data[0]
                        except EOFError:
                            break

                        for band in range(num_bands):
                            occup = 1.0/self.dense_usage_count[densept] * dense_occupation_dict[densept][band]
                            occups_dict[tuple(coarsept)][band]+=occup
                            total_weighted_occupation += occup
                        
            #import pdb; pdb.set_trace()
            
            print(40*"-")
        
            print(f"\nWriting occupations to {occupationfilename} . . .\n")
        
            print(40*"-")
            
            writer = csv.writer(occupationfile,delimiter=",")
            
            if band_labels:
                writer.writerow(["k1","k2","k3",*band_labels])
            
            else:
                labels = [f"f_{band}" for band in range(0,num_bands)]
                writer.writerow(["k1","k2","k3",*labels])
            
            #write occupations to file
            for coarsept,occups in tqdm(occups_dict.items()):
                writer.writerow([*coarsept,*occups])

            #Ensure normalization
            assert abs(total_weighted_occupation - 1.0) < 1e-6, f"Normalization of the occupations is not correct! Got {total_weighted_occupation} instead of 1.0!"

            print(f"\n Normalization is correct! Got: {total_weighted_occupation}")

            print("\n DONE! \n")
    
    
    def plot_grids(self,neighbourfilename,coarse=True,dense=True):
        
        """
        Tool to visualize the different grids. Doesn't work well for large grids, due to many points. 
        
        Can possible add an option to view a portion of it, just like an isosurface. 
        
        
        """
        
        neighbourfile = open(neighbourfilename,"r")     
        reader = csv.reader(neighbourfile,delimiter=",")
        next(reader) # skip header
        
        #check if all kpoints in each direction are available, and if not, calculate
        

        # for the future
        
        #get all unique k-points for coarse grid and/or dense grid, if true
        
        coarse_different = set()
        dense_different = set()
        
        for i,ln in enumerate(reader):
                
            ln = np.genfromtxt(ln)
            
            coarsept,densept = tuple(ln[0:3]),tuple(ln[3:6])
            
            if coarse is True:
                coarse_different.add(coarsept)
            
            if dense is True:
                dense_different.add(densept)
        
        #collect points in three different sublists
        
        coarse_arrays = defaultdict(list)
        for k in coarse_different:
            for i in range(0,3):
                coarse_arrays[i].append(k[i])

                
        dense_arrays = defaultdict(list)
        for k in dense_different:
            for i in range(0,3):
                dense_arrays[i].append(k[i])
        

        #plot fig 
        
        fig = plt.figure(figsize = (8,8))
        ax = plt.axes(projection='3d')
        ax.grid()
        ax.scatter(dense_arrays[0],dense_arrays[1],dense_arrays[2],s=15,color="green",alpha=0.5,label=f"dense grid") #{len(dense_arrays[0])}x{len(dense_arrays[1])}x{len(dense_arrays[2])}",alpha=0.5) --> need to implement a way to get all different points 
        ax.scatter(coarse_arrays[0],coarse_arrays[1],coarse_arrays[2],s=25,color="red",label=f"coarse grid") #{len(coarse_arrays[0])}x{len(coarse_arrays[1])}x{len(coarse_arrays[2])}",color="red")
        
        ax.legend() 
        
        plt.show()


def grid_match_test():
    
    """
    
    Test the functionalities with dummy grids. 
    
    """
    
    #test with dummy grids
    yambogrid = np.linspace(0.0,1.0,4,endpoint=False)
    yambogrid3D = np.array(list(product(yambogrid,yambogrid,yambogrid)))
    pertgrid = np.linspace(0.0,1.0,12,endpoint=False)
    pertgrid3D = np.array(list(product(pertgrid,pertgrid,pertgrid)))

    #generate dummy occupation dictionary
    pert_occupations = np.random.rand(6,len(pertgrid3D))
    pert_occupations /= pert_occupations.sum() #normalize
    
    pert_occupations_dict = {}
    
    for pertpt,pertf in zip(pertgrid3D,pert_occupations.T[::-1]):
        pert_occupations_dict[tuple(pertpt.tolist())] = pertf
        
    #retrieve nearest neighbours for the yambo grid
    grid = GridMatch()
    grid.get_nearest_neighbours(yambogrid3D,pertgrid3D,filename="neighbours_occups_test_160625.csv")
    
    #match occupations
    grid.gen_coarse_occupations(pert_occupations_dict,neighbourfilename="neighbours_occups_test_160625.csv")
    
    #plot lattices 
    grid.plot_grids("neighbours_occups_test_160625.csv")

if __name__ == "__main__":
    
    grid_match_test()
    