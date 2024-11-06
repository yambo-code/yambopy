import numpy as np
from scipy.fft import ifftn,fftn,fft2

class FourierTransformer:
    """
    Class to handle FFT transformation. 
    """    
    def __init__(self, dlat,rlat,qpoints,fft_grid_size): #fft_grid_size can 3D,2D array
        '''
        yambopy FFT class. Perform the np.fft transform of 2D and 3D grids
        fft_grid_size can be a tuple, array or numpy array.
        '''
        
        #mode can be either real to direct or direct2real 'd2r'
        self.dlat = dlat # direct lattice
        self.rlat = rlat #reciprocal lattice
        self.qpoints = np.array(qpoints)
        self.fft_grid_size = fft_grid_size
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
        x,y,z = np.linspace(0,self.fft_grid_size[0],1),np.linspace(0,self.fft_grid_size[1],1),np.linspace(0,self.fft_grid_size[2],1)
        if self.dim == 2:
            X, Y = np.meshgrid(x, y)
            return X, Y
            #space_grid = np.column_stack((x.ravel(), y.ravel()))
        elif self.dim == 3:
            X, Y, Z = np.meshgrid(x, y, z)
            #space_grid = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
            return X, Y, Z
    def real_to_reciprocal(self, data_real):
        # Perform the Fourier transform to get the function in reciprocal space
        data_reciprocal = fftn(data_real, s = self.fft_grid_size)
        return data_reciprocal