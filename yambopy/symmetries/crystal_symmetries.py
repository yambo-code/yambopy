import numpy as np
import spgrep
import spglib
from yambopy.dbs.latticedb import YamboLatticeDB

class Crystal_Symmetries:
    def __init__(self, latdb, magnetic_moments=None, tol=1e-5):
        """
        Class to handle  crystal symmetries
        #
        Given a lattice database, finds all symmetries of a  crystal.
        Works irrespective of whether the SAVE has symmetries or not.
        #
        Parameters:
        -----------
        latdb : object
            Lattice database object with required attributes:
        magnetic_moments : array (natom) or (natom,3) or None, optional
            magnetic_moments of each atom (default: None)
            In collinear case, natom is sufficient, but in non-collinear 
            case provide (mx, my, mz) for each atom.
        tol : float, optional
            Symmetry tolerance (default: 1e-5)
        Members:
        --------
        dataset : dict
            The full symmetry dataset from spglib
        spacegroup_type : dict
            Space group type information
        rotations : ndarray
            Rotation matrices in Cartesian coordinates (n_sym, 3, 3)
        translations : ndarray
            Translation vectors in Cartesian coordinates (n_sym, 3)
        international_symbol : str
            International symbol of the space group
        hall_symbol : str
            Hall symbol of the space group
        wyckoffs : list
            Wyckoff letters for each atom
        pointgroup : str
            Point group symbol
        spacegroup_number : int
            International space group number
        """
        # Get lattice and atomic information
        lattice = latdb.lat
        numbers = np.rint(latdb.atomic_numbers).astype(int)
        lat_inv = np.linalg.inv(lattice)
        positions = latdb.car_atomic_positions @ lat_inv
        positions = positions - np.floor(positions)
        #
        self.magnetic = not (magnetic_moments is None)
        #
        # Get symmetry dataset
        if self.magnetic :
            cell = (lattice, positions, numbers, magnetic_moments)
            print("Currently Magnetic systems not supported.")
            exit()
        else :
            cell = (lattice, positions, numbers)
            self.dataset = spglib.get_symmetry_dataset(cell, symprec=tol)
            self.spacegroup_type = spglib.get_spacegroup_type(self.dataset.hall_number)
        #
        # Convert rotation and fractional translations to Cartesian units
        self.rotations = lattice.T[None,:,:] @ self.dataset.rotations @ lat_inv.T[None,:,:]
        self.translations = self.dataset.translations @ lattice
        #
        self.international_symbol = self.dataset.international
        self.hall_symbol = self.dataset.hall
        self.wyckoffs = self.dataset.wyckoffs
        self.pointgroup = self.dataset.pointgroup
        self.pointgroup_schoenflies = self.spacegroup_type.pointgroup_schoenflies
        self.spacegroup_number = self.dataset.number


## Test
if __name__ == "__main__":
    import os
    np.set_printoptions(suppress=True)
    latdb = YamboLatticeDB.from_db_file(os.path.join('.', 'SAVE'))
    symm = Crystal_Symmetries(latdb,tol=1e-4)
    print(symm.wyckoffs)
    print(symm.pointgroup_schoenflies)

