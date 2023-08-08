import numpy as np
from scipy.fft import ifftn,fftn

class FourierTransformer:
    """
    Class to handle FFT transformation. 
    """    
    def __init__(self, dlat,rlat,qpoints,fft_grid_size, mode='r2d'): #fft_grid_size can 3D,2D array
        #mode can be either real to direct or direct2real 'd2r'
        self.dlat = dlat # direct lattice
        self.rlat = rlat #reciprocal lattice
        self.qpoints = np.array(qpoints)
        self.fft_grid_size = np.array(fft_grid_size)
        self.mode = mode
        # Compute the number of G-vectors based on the fft_grid_size
        self.num_g_vectors = np.prod(fft_grid_size)

        # Check the dimensionality of qpoints and adjust for 2D case if necessary
        if self.qpoints.shape[1] == 2:
            self.dim = 2
        elif self.qpoints.shape[1] == 3:
            self.dim = 3
        # elif self.qpoints.shape[1] == 1:
        #     self.dim = 1            
        else:
            raise ValueError("Invalid dimensionality for qpoints. Must be 2D or 3D.")        

        if self.dim != len(fft_grid_size):
            raise ValueError("Mismatch between dimensionality of qpoints and fft_grid_size.")
        #number of G-vecotrs based on the qpoints

    def get_grid(self):
        # Create the grid for real space coordinates
        if self.dim == 2:
            x, y = np.meshgrid(np.arange(self.fft_grid_size[0]+1), np.arange(self.fft_grid_size[1]+1))
            space_grid = np.column_stack((x.ravel(), y.ravel()))
        elif self.dim == 3:
            x, y, z = np.meshgrid(np.arange(self.fft_grid_size[0]+1), np.arange(self.fft_grid_size[1]+1), np.arange(self.fft_grid_size[2]+1))
            space_grid = np.column_stack((x.ravel(), y.ravel(), z.ravel()))

        return space_grid        

    def get_rgrid(self):
        # Create the grid with points in units of the reciprocal lattice vectors
        if self.dim == 2:
            x, y = np.mgrid[0:1:self.fft_grid_size[0]*1j, 0:1:self.fft_grid_size[1]*1j]
            space_grid = np.column_stack((x.ravel(), y.ravel()))
        elif self.dim == 3:
            x, y, z = np.mgrid[0:1:self.fft_grid_size[0]*1j, 0:1:self.fft_grid_size[1]*1j,0:1:self.fft_grid_size[2]*1j]
            space_grid = np.column_stack((x.ravel(), y.ravel(), z.ravel()))        
        rspace_grid = np.dot(space_grid,self.rlat)

        return rspace_grid

    def get_dgrid(self):
        # Create the grid with points in units of the reciprocal lattice vectors
        if self.dim == 2:
            x, y = np.mgrid[0:1:self.fft_grid_size[0]*1j, 0:1:self.fft_grid_size[1]*1j]
            space_grid = np.column_stack((x.ravel(), y.ravel()))
        elif self.dim == 3:
            x, y, z = np.mgrid[0:1:self.fft_grid_size[0]*1j, 0:1:self.fft_grid_size[1]*1j,0:1:self.fft_grid_size[2]*1j]
            space_grid = np.column_stack((x.ravel(), y.ravel(), z.ravel()))        
        dspace_grid = np.dot(space_grid,self.lat)

        return dspace_grid

    def find_nearest_grid_point(self, qpoint, real_space_grid):
        # Find the nearest grid point in the real space grid for a given q-point
        rounded_q_point = np.around(qpoint).astype(int)
        nearest_grid_point_idx = np.argmin(np.linalg.norm(real_space_grid - rounded_q_point, axis=1))
        return nearest_grid_point_idx,real_space_grid[nearest_grid_point_idx]

    def transform_to_real_space(self, function_values):
        # Perform the inverse Fourier transform to get the function in real space
        real_space_data = ifftn(function_values).real
        return real_space_data

    def transform_to_reciprocal_space(self, real_space_function):
        # Perform the Fourier transform to get the function in reciprocal space
        reciprocal_space_data = fftn(real_space_function)
        return reciprocal_space_data

    def q_point_to_real_space(self, qpoint):
        # Convert a q-point to real space coordinates in the unit cell
        real_space_point = np.dot(qpoint, self.rlat)
        return real_space_point

    def q_points_to_real_space(self):
        # Convert all q-points to real space coordinates in the unit cell
        real_space_points = np.dot(self.qpoints, self.rlat)
        return real_space_points

    def real_space_to_q_point(self, real_space_point):
        # Convert a real space point in the unit cell to a q-point
        qpoint = np.linalg.solve(self.rlat_lattice.T, real_space_point)
        return qpoint

    def real_space_to_q_points(self, real_space_points):
        # Convert all real space points in the unit cell to q-points
        qpoints = np.linalg.solve(self.rlat_lattice.T, real_space_points.T).T
        return qpoints

    def perform_fourier_transform(self, function_values):
        if self.mode == "r2d":
            # Transform function from reciprocal to real space
            real_space_functions = []
            space_grid = self.get_grid()
            rspace_grid = self.get_rgrid()
            print(rspace_grid)
            for qpoint in self.qpoints:
                real_space_point = self.q_point_to_real_space(qpoint)
                print(real_space_point)
                nearest_grid_point_idx, nearest_grid_point = self.find_nearest_grid_point(real_space_point, rspace_grid)
                print(nearest_grid_point[0],nearest_grid_point[1],nearest_grid_point[2])
                function_values_at_nearest_point = function_values[nearest_grid_point[0], nearest_grid_point[1], nearest_grid_point[2]]
                print(function_values_at_nearest_point)
                real_space_function = self.transform_to_real_space(function_values_at_nearest_point)
                real_space_functions.append(real_space_function)
            return real_space_functions
        elif self.transform_direction == "d2r":
            # Transform function from real to reciprocal space
            reciprocal_space_functions = []
            for real_space_function in function_values:
                reciprocal_space_function = self.transform_to_reciprocal_space(real_space_function)
                reciprocal_space_functions.append(reciprocal_space_function)
            return reciprocal_space_functions
        else:
            raise ValueError("Invalid transform_direction. Must be 'reciprocal_to_real' or 'real_to_reciprocal'.")

