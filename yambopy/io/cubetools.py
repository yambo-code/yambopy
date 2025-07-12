import numpy as np

def write_cube(filename, data, lat_vec, atom_pos, atomic_num, origin=np.zeros(3), header=''):
    """
    Generate a Gaussian cube for a real data (Nx, Ny, Nz)
    filename : Name of the file
    data : real data with dimension (Nx, Ny, Nz)
    lat_vec : lattice vectors ai = a[:,i]
    origin : origin of the cell 
    atom_pos: atomic positions in cart units (bohr)
    atomic_num: atomic numbers of elements
    """
    ## atomic positions must be in [0,1)
    apos_crys = atom_pos@lat_vec.T
    apos_crys = apos_crys-np.floor(apos_crys)
    apos_crys = (apos_crys+1e-6)%1
    apos_new = apos_crys@np.linalg.inv(lat_vec.T)

    with open(filename, "w") as cube:
        cube.write("# Cube file generated using YamboPy\n")
        cube.write('# '+header.strip()+'\n')
        cube.write("%d     %.6f     %.6f      %.6f\n" \
                   %(len(atom_pos), origin[0], origin[1],origin[2] ))
        for i in range(3):
            Ntmp = data.shape[i]
            cube.write("%d     %.6f     %.6f      %.6f\n" \
                   %(Ntmp, lat_vec[0,i]/Ntmp, lat_vec[1,i]/Ntmp,lat_vec[2,i]/Ntmp))

        for i in range(len(atom_pos)):
            cube.write("%d     %d     %.6f     %.6f      %.6f\n" \
                   %(atomic_num[i], 0.0, apos_new[i,0], apos_new[i,1], apos_new[i,2]))

        for ix in range(data.shape[0]):
            for iy in range(data.shape[1]):
                for iz in range(data.shape[2]):
                    cube.write("%.6f      " %(data[ix,iy,iz]))
                    if iz % 6 == 5 : cube.write("\n")

