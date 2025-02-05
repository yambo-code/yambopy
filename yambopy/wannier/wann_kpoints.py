from yambopy.lattice import red_car
import numpy as np
from scipy.spatial import cKDTree
class KPointGenerator():
    def __init__(self):
        self.k = None
        self.nkpoints = None
        self.weights = None
        self.red_kpoints = None
        self.car_kpoints = None
        self.k_tree = None

    def generate(self):
        """Abstract method to generate k-points."""
        raise NotImplementedError("This method must be implemented in subclasses.")

    def validate(self):
        """Validate the generated k-points."""
        if self.k is None or self.weights is None:
            raise ValueError("k-points or weights have not been generated.")
        if self.red_kpoints or self.car_kpoints is None:
            raise ValueError("red and/or cartesian k-points not initialized.")
    def export(self, filename):
        """Export k-points to a file."""
        np.savetxt(filename, self.k, header="k-points")

    def fold_into_bz_Gs(self, k_points, bz_range=(-0.5, 0.5), reciprocal_vectors=None):
        """
        Fold k-points into the first Brillouin Zone and determine the reciprocal lattice vectors G needed.

        Parameters:
        - k_points: Array of k-points in k-space (shape: (..., 3)).
        - bz_range: Tuple indicating the range of the BZ, default is (-0.5, 0.5) for each direction.
        - reciprocal_vectors: A list or matrix of reciprocal lattice vectors defining the BZ.

        Returns:
        - folded_k_points: The folded k-points within the BZ (shape: same as k_points).
        - G_vectors: The reciprocal lattice vectors that fold the k-points into the BZ (shape: same as k_points).
        """
        k_points = np.array(k_points)  # Ensure input is an array

        # Determine the G-vector multipliers for folding
        G_multiplier = np.floor((k_points - bz_range[0]) / (bz_range[1] - bz_range[0]))
        
        # Calculate the G_vectors
        if reciprocal_vectors is not None:
            reciprocal_vectors = np.array(reciprocal_vectors)  # Ensure reciprocal_vectors is an array
            G_vectors = np.tensordot(G_multiplier, reciprocal_vectors, axes=([1], [0]))
        else:
            # Assume a cubic lattice with unit cell length of 1
            G_vectors = G_multiplier * (bz_range[1] - bz_range[0])

        # Fold the k-points into the BZ
        folded_k_points = k_points - G_vectors
        # Correct for points exactly on the upper bound of the BZ
        mask_upper_bound = (folded_k_points == bz_range[0]) & (G_vectors >= 1.0)
        folded_k_points[mask_upper_bound] += (bz_range[1] - bz_range[0])
        G_vectors[mask_upper_bound] -= (bz_range[1] - bz_range[0])
        # Negate G_vectors
        G_vectors = -G_vectors
        # Handle points at the upper bound of the BZ
        mask_at_upper_bound = (folded_k_points == bz_range[1])
        folded_k_points[mask_at_upper_bound] -= (bz_range[1] - bz_range[0])
        G_vectors[mask_at_upper_bound] = 0

        return folded_k_points, G_vectors


    def find_closest_kpoint(self, points):

        # Convert points to a numpy array
        points = np.atleast_2d(points)  # Ensure points is a 2D array (N, 3)
        
        # Calculate distances using broadcasting and periodic boundary conditions
        #delta = (self.k[np.newaxis, :, :] - points[:, np.newaxis, :] + 0.5) % 1 - 0.5
        distance, closest_indices = self.k_tree.query(points, k=1)
        #distances = np.linalg.norm(delta, axis=2)  # Shape (N, nkpoints)
        
        # Find the index of the minimum distance for each point
        #closest_indices = np.argmin(distances, axis=1)  # Shape (N,)
        
        # Return the appropriate type
        return closest_indices if points.shape[0] > 1 else int(closest_indices[0])


    def get_kq_tables(self,qmpgrid, sign = '+'):
        kplusq_table = np.zeros((self.nkpoints,qmpgrid.nkpoints),dtype=int)
        kminusq_table = np.zeros((self.nkpoints,qmpgrid.nkpoints), dtype=int)
        # Assign to tables

        _,kplusq_table = self.get_kpq_grid(qmpgrid)
        _,kminusq_table = self.get_kmq_grid(qmpgrid,sign)

        return kplusq_table, kminusq_table
    

    def get_kindices_fromq(self,qmpgrid):
        q_points = qmpgrid.k  # Extract q-points array

        # Apply the function self.find_closest_kpoint to all q-points in a vectorized manner
        kindices_fromq = np.apply_along_axis(self.find_closest_kpoint, -1, q_points)

        # Assign the result directly
        self.kindices_fromq = kindices_fromq

        return kindices_fromq

    def fold_into_bz(self,points):
        'Fold a point in the first BZ defined in the range [-0.5,0.5]'
        # Applying the modulo operation to shift points within the range [-0.5, 0.5]
        folded_points = np.mod(points + 0.5, 1.0) - 0.5
        # Correcting points where original points were exactly 0.5 to remain 0.5
        folded_points[(points == 0.5)] = 0.5
        return folded_points
    
    def get_plaquette(self, nx, ny, nz, dir = 2):
        if (dir == 2):
            dir1 = np.array([1,0,0])
            dq1 = 1/nx
            dir2 = np.array([0,1,0])
            dq2 = 1/ny
            nps = nx*ny
        if (dir == 0):
            dir1 = np.array([0,1,0])
            dq1 = 1/ny
            dir2 = np.array([0,0,1])
            dq2 = 1/nz
            nps = ny*nz
        if (dir == 1):
            dir1 = np.array([1,0,0])
            dq1 = 1/nx
            dir2 = np.array([0,0,1])            
            dq2 = 1/nz
            nps = nx*nz
        counter = 0
        #here I need to use the q grid and then apply -b/2
        qplaquette_grid = np.zeros((nps, 4), dtype=int)
        for iq, q in enumerate(self.k):     # needs improving with array casting
            if (q[dir]==0.0):
                tmp_qp1, tmp_qp1Gvec = self.fold_into_bz_Gs(q + dq1*dir1)
                tmp_qp1p2, tmp_qp1p2Gvec = self.fold_into_bz_Gs(q + dq1*dir1+ dq2*dir2)
                tmp_qp2, tmp_qp2Gvec = self.fold_into_bz_Gs(q + dq2*dir2)  
                idxqp1 = self.find_closest_kpoint(tmp_qp1)
                idxqp1p2 = self.find_closest_kpoint(tmp_qp1p2)
                idxqp2 = self.find_closest_kpoint(tmp_qp2)
                qplaquette_grid[counter] = [iq, idxqp1, idxqp1p2, idxqp2]
                counter +=1            

        self.qplaquette_grid = qplaquette_grid