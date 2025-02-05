from yambopy.lattice import red_car
from yambopy.units import bohr2ang
import numpy as np
from yambopy.wannier.wann_io import NNKP
from yambopy.wannier.wann_kpoints import KPointGenerator
from yambopy.units import ang2bohr
import numpy as np
from scipy.spatial import cKDTree

class tb_Monkhorst_Pack(KPointGenerator):
    def __init__(self, grid_shape,latdb, shift=np.array([0.0,0.0,0.0])):
        super().__init__()
        self.grid_shape = grid_shape
        self.latdb = latdb
        self.rlat = self.latdb.rlat*2*np.pi*ang2bohr
        self.shift = shift
    def generate(self):
        # Create grid points for n1, n2, and n3
        NGX, NGY, NGZ = self.grid_shape        
        n1 = np.arange(NGX)
        n2 = np.arange(NGY)
        n3 = np.arange(NGZ)

        e1 = np.array([1.0, 0.0, 0.0])
        e2 = np.array([0.0, 1.0, 0.0])  
        e3 = np.array([0.0, 0.0, 1.0])  
        # Create a meshgrid for all combinations of n1, n2, and n3
        n1, n2, n3 = np.meshgrid(n1, n2, n3, indexing='ij')
        b1,b2,b3 = self.rlat
        # Calculate the k-points
        k_points = (n1[:, :, :, np.newaxis] / NGX) * b1 + (n2[:, :, :, np.newaxis] / NGY) * b2 + (n3[:, :, :, np.newaxis] / NGZ) * b3 + self.rlat@self.shift
        red_points = (n1[:, :, :, np.newaxis] / NGX) * e1 + (n2[:, :, :, np.newaxis] / NGY) * e2 + (n3[:, :, :, np.newaxis] / NGZ) * e3 + self.shift
        # Flatten the grid into a list of k-points
        self.k = red_points.reshape(-1, 3)
        self.red_kpoints = self.k
        self.car_kpoints = k_points.reshape(-1,3)
        self.nkpoints = len(self.k)
        self.k_tree = cKDTree(self.k)        

    def get_kmq_grid(self,qmpgrid, sign):
        # if not isinstance(qmpgrid, NNKP_Grids):
        #     raise TypeError('Argument must be an instance of NNKP_Grids')
        #here I need to use the k-q grid and then apply -b/2
        # Prepare dimensions
        if sign not in ["+", "-"]:
            raise ValueError("Invalid sign option. Choose either '+' for a-b or '-' for b-a")
        
        nkpoints = self.nkpoints
        nqpoints = qmpgrid.nkpoints
     
        # Broadcast k and q grids to shape (nkpoints, nqpoints, 3)
        k_grid = np.expand_dims(self.k, axis=1)  # Shape (nkpoints, 1, 3)
        q_grid = np.expand_dims(qmpgrid.k, axis=0)  # Shape (1, nqpoints, 3)
        if sign == "+":
            kq_diff = k_grid - q_grid  # Shape (nkpoints, nqpoints, 3)
        elif sign == "-":
            kq_diff = -k_grid + q_grid  # Shape (nkpoints, nqpoints, 3)

        # Fold into the Brillouin Zone and get G-vectors
        kmq_folded, Gvec = self.fold_into_bz_Gs(kq_diff.reshape(-1, 3),bz_range=(0.0,1.0))  # Flatten for batch processing
        kmq_folded = kmq_folded.reshape(nkpoints, nqpoints, 3)
        Gvec = Gvec.reshape(nkpoints, nqpoints, 3)

        # Find closest k-points for all points
        closest_indices = self.find_closest_kpoint(kmq_folded.reshape(-1, 3)).reshape(nkpoints, nqpoints)
        # Populate the grids
        kmq_grid = kmq_folded  # Shape (nkpoints, nqpoints, nnkpts, 3)
        kmq_grid_table = np.stack(
            [
                np.arange(nkpoints)[:, None].repeat(nqpoints, axis=1),  # ik
                closest_indices,  # idkmq
                Gvec[..., 0].astype(int),  # Gx
                Gvec[..., 1].astype(int),  # Gy
                Gvec[..., 2].astype(int)   # Gz
            ],
            axis=-1
        ).astype(int)  # Shape (nkpoints, nqpoints, 5)
        self.kmq_grid = kmq_grid
        self.kmq_grid_table = kmq_grid_table        
        
        return kmq_grid, kmq_grid_table

    def get_kmqpdq_grid(self,qmpgrid):

        nkpoints = self.nkpoints
        nqpoints = qmpgrid.nkpoints

        # Broadcast k and q grids to shape (nkpoints, nqpoints, 3)
        k_grid = np.expand_dims(self.k, axis=1)  # Shape (nkpoints, 1, 3)
        q_grid = np.expand_dims(qmpgrid.k, axis=0)  # Shape (1, nqpoints, 3)
        NX, NY, NZ = self.grid_shape
        spacing = np.array([[1/NX, 0.0, 0.0],[0.0, 1/NY, 0.0],[0.0, 0.0, 1/NZ]], dtype=np.float128)
        kmqdq_grid = []
        kmqdq_grid_table = []
        for idq, dq in enumerate(spacing):
            kq_diff = k_grid - q_grid+dq  # Shape (nkpoints, nqpoints, 3)

            # Fold into the Brillouin Zone and get G-vectors
            kmq_folded, Gvec = self.fold_into_bz_Gs(kq_diff.reshape(-1, 3),bz_range=(0.0,1.0))  # Flatten for batch processing
            kmq_folded = kmq_folded.reshape(nkpoints, nqpoints, 3)
            Gvec = Gvec.reshape(nkpoints, nqpoints, 3)

            # Find closest k-points for all points
            closest_indices = self.find_closest_kpoint(kmq_folded.reshape(-1, 3)).reshape(nkpoints, nqpoints)
            # Populate the grids
            kmq_grid = kmq_folded  # Shape (nkpoints, nqpoints, nnkpts, 3)
            kmq_grid_table = np.stack(
                [
                    np.arange(nkpoints)[:, None].repeat(nqpoints, axis=1),  # ik
                    closest_indices,  # idkmq
                    Gvec[..., 0].astype(int),  # Gx
                    Gvec[..., 1].astype(int),  # Gy
                    Gvec[..., 2].astype(int)   # Gz
                ],
                axis=-1
            ).astype(int)  # Shape (nkpoints, nqpoints, 5)
            kmqdq_grid.append(kmq_grid)
            kmqdq_grid_table.append(kmq_grid_table)        
        
        self.kmqdq_grid = kmqdq_grid
        self.kmqdq_grid_table = kmqdq_grid_table

        return kmqdq_grid, kmqdq_grid_table
        
    def get_kpq_grid(self, qmpgrid):

        nkpoints = self.nkpoints
        nqpoints = qmpgrid.nkpoints
        
        # Broadcast k and q grids to shape (nkpoints, nqpoints, 3)
        k_grid = np.expand_dims(self.k, axis=1)  # Shape (nkpoints, 1, 3)
        q_grid = np.expand_dims(qmpgrid.k, axis=0)  # Shape (1, nqpoints, 3)
        kq_add = k_grid + q_grid  # Shape (nkpoints, nqpoints, 3)

        # Fold into the Brillouin Zone and get G-vectors
        kpq_folded, Gvec = self.fold_into_bz_Gs(kq_add.reshape(-1, 3),bz_range=(0.0,1.0))  # Flatten for batch processing
        kpq_folded = kpq_folded.reshape(nkpoints, nqpoints, 3)
        Gvec = Gvec.reshape(nkpoints, nqpoints, 3)

        # Find closest k-points for all points
        closest_indices = self.find_closest_kpoint(kpq_folded.reshape(-1, 3)).reshape(nkpoints, nqpoints)
        # Populate the grids
        kpq_grid = kpq_folded  # Shape (nkpoints, nqpoints, nnkpts, 3)
        kpq_grid_table = np.stack(
            [
                np.arange(nkpoints)[:, None].repeat(nqpoints, axis=1),  # ik
                closest_indices,  # idxkp
                Gvec[..., 0].astype(int),  # Gx
                Gvec[..., 1].astype(int),  # Gy
                Gvec[..., 2].astype(int)   # Gz
            ],
            axis=-1
        ).astype(int)  # Shape (nkpoints, nqpoints, 5)

        self.kpq_grid = kpq_grid
        self.kpq_grid_table = kpq_grid_table

        return kpq_grid, kpq_grid_table

    def get_qpb_grid(self, qmpgrid: 'tb_Monkhorst_Pack'):
        '''
        For each q belonging to the Qgrid return Q+B and a table with indices
        containing the q index the q+b folded into the BZ and the G-vectors
        '''
        if not isinstance(qmpgrid, tb_Monkhorst_Pack):
            raise TypeError('Argument must be an instance of tb_Monkhorst_Pack')
        
        # here I should work only with qmpgrid and its qmpgrid.b_grid
        qpb_grid = np.zeros((qmpgrid.nkpoints, qmpgrid.nnkpts, 3))
        qpb_grid_table = np.zeros((qmpgrid.nkpoints, qmpgrid.nnkpts, 5), dtype = int)
        for iq, q in enumerate(qmpgrid.k):
            for ib, b in enumerate(qmpgrid.b_grid[qmpgrid.nnkpts*iq:qmpgrid.nnkpts*(iq+1)]):
                tmp_qpb, tmp_Gvec = qmpgrid.fold_into_bz_Gs(q+b)
                idxqpb = self.find_closest_kpoint(tmp_qpb)
                qpb_grid[iq, ib] = tmp_qpb
                # here it should be tmp_Gvec, but with yambo grid I have inconsistencies because points are at 0.75
                qpb_grid_table[iq,ib] = [iq, idxqpb, int(qmpgrid.iG[ib+qmpgrid.nnkpts*iq,0]), int(qmpgrid.iG[ib+qmpgrid.nnkpts*iq,1]), int(qmpgrid.iG[ib+qmpgrid.nnkpts*iq,2])]
        
        self.qpb_grid = qpb_grid
        self.qpb_grid_table = qpb_grid_table

    def get_kpbover2_grid(self, qmpgrid: 'tb_Monkhorst_Pack'):
        if not isinstance(qmpgrid, tb_Monkhorst_Pack):
            raise TypeError('Argument must be an instance of tb_Monkhorst_Pack')

        # Reshape b_grid for vectorized addition
        b_grid = self.b_grid.reshape(self.nkpoints, self.nnkpts, 3)  # Shape (nkpoints, nnkpts, 3)

        # Add k and b_grid with broadcasting
        k_expanded = self.k[:, np.newaxis, :]  # Shape (nkpoints, 1, 3)
        combined_kb = k_expanded + b_grid  # Shape (nkpoints, nnkpts, 3)
        # Fold into the BZ for all k + b combinations
        folded_kb, Gvec = self.fold_into_bz_Gs(combined_kb.reshape(-1, 3),bz_range=(0.0,1.0))  # Flatten first two dims
        folded_kb = folded_kb.reshape(self.nkpoints, self.nnkpts, 3)
        Gvec = Gvec.reshape(self.nkpoints, self.nnkpts, 3)
        # Find closest kpoints
        idxkpbover2 = self.find_closest_kpoint(folded_kb.reshape(-1, 3)).reshape(self.nkpoints, self.nnkpts)
        # # Construct results

        self.kpbover2_grid = folded_kb
        self.kpbover2_grid_table = np.stack([
            np.repeat(np.arange(self.nkpoints)[:, np.newaxis], self.nnkpts, axis=1),  # ik
            idxkpbover2,  # Closest kpoint indices
            Gvec[..., 0].astype(int),  # Gx
            Gvec[..., 1].astype(int),  # Gy
            Gvec[..., 2].astype(int)   # Gz
        ], axis=-1)  # Shape (nkpoints, nnkpts, 5)


    def get_kmqmbover2_grid(self, qmpgrid: 'tb_Monkhorst_Pack'):       # need to improve this one
        if not isinstance(qmpgrid, tb_Monkhorst_Pack):
            raise TypeError('Argument must be an instance of tb_Monkhorst_Pack')

        # Prepare dimensions
        nkpoints = self.nkpoints
        nqpoints = qmpgrid.nkpoints
        nnkpts = qmpgrid.nnkpts

        # Broadcast k and q grids to shape (nkpoints, nqpoints, 3)
        k_grid = np.expand_dims(self.k, axis=1)  # Shape (nkpoints, 1, 3)
        q_grid = np.expand_dims(qmpgrid.k, axis=0)  # Shape (1, nqpoints, 3)
        kq_diff = k_grid - q_grid  # Shape (nkpoints, nqpoints, 3)

        # Reshape b_grid to be specific to each k-point
        b_grid = self.b_grid.reshape(nkpoints, nnkpts, 3)  # Shape (nkpoints, nnkpts, 3)

        # Calculate k - q - b for all combinations
        kqmbover2 = (
            kq_diff[:, :, np.newaxis, :]  # Shape (nkpoints, nqpoints, 1, 3)
            - b_grid[:, np.newaxis, :, :]  # Shape (nkpoints, 1, nnkpts, 3)
        )  # Final Shape (nkpoints, nqpoints, nnkpts, 3)

        # Fold into the Brillouin Zone and get G-vectors
        kqmbover2_folded, Gvec = self.fold_into_bz_Gs(kqmbover2.reshape(-1, 3),bz_range=(0.0,1.0))  # Flatten for batch processing
        kqmbover2_folded = kqmbover2_folded.reshape(nkpoints, nqpoints, nnkpts, 3)
        Gvec = Gvec.reshape(nkpoints, nqpoints, nnkpts, 3)

        # Find closest k-points for all points
        closest_indices = self.find_closest_kpoint(kqmbover2_folded.reshape(-1, 3)).reshape(nkpoints, nqpoints, nnkpts)
        # Populate the grids
        self.kmqmbover2_grid = kqmbover2_folded  # Shape (nkpoints, nqpoints, nnkpts, 3)
        self.kmqmbover2_grid_table = np.stack(
            [
                np.arange(nkpoints)[:, None, None].repeat(nqpoints, axis=1).repeat(nnkpts, axis=2),  # ik
                closest_indices,  # idxkmqmbover2
                Gvec[..., 0].astype(int),  # Gx
                Gvec[..., 1].astype(int),  # Gy
                Gvec[..., 2].astype(int)   # Gz
            ],
            axis=-1
        ).astype(int)  # Shape (nkpoints, nqpoints, nnkpts, 5)