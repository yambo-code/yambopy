import numpy as np
from qepy import PwIn
from itertools import product
import copy
import math

pw = PwIn('scf/bn.scf')
latvec = pw.cell_parameters
basis  = int(pw.system['nat'])
atoms  = pw.atoms
R = [4,1,1]

# [1] Generate atomic positions in the supercell
#===============================================================================
def supercell(atoms,latvec,basis,R):
    sup_size = R[0]*R[1]*R[2]
    latvec = np.array(latvec)
    new_latvec = np.array([R[0]*latvec[0],R[1]*latvec[1],R[2]*latvec[2]])
    atoms = np.array([atom[1] for atom in atoms])
    #new_atoms[cell][basis][direction]
    new_atoms = np.array([atoms for n in range(sup_size)])
    for nz,ny,nx,b in product(range(R[2]),range(R[1]),range(R[0]),range(basis)):
        cell=nx+ny*R[0]+nz*R[0]*R[1]
        new_atoms[cell,b,0]=(new_atoms[cell,b,0]+nx)/R[0]
        new_atoms[cell,b,1]=(new_atoms[cell,b,1]+ny)/R[1]
        new_atoms[cell,b,2]=(new_atoms[cell,b,2]+nz)/R[2]
    #new_atoms[super_basis][directions]
    new_atoms=new_atoms.reshape(basis*sup_size,3)
    return new_latvec, new_atoms
#
# [2] Generate a list of [ atom_type [position] ] for the supercell input file
#===============================================================================
def atoms_input(uc_input, uc_basis, sc_positions, sc_size):
    sc_basis=sup_size*uc_basis
    positions_input = sc_positions.tolist()
    elements_input  = [[uc_input[i][0] for i in range(uc_basis)] for j in range(sc_size)]
    elements_input  = [ item for sublist in elements_input for item in sublist ]
    atoms_input     = [[elements_input[i], positions_input[i]] for i in range(sc_basis)]
    return atoms_input
#===============================================================================

s_latvec,s_atoms= supercell(atoms,latvec,basis,R)
print s_atoms
sup_size = R[0]*R[1]*R[2]
new_kpoints = [math.ceil(pw.kpoints[0]/R[0]), math.ceil(pw.kpoints[1]/R[1]), math.ceil(pw.kpoints[2]/R[2])]
pw_s = copy.deepcopy(pw)
pw_s.atoms = atoms_input(pw.atoms,basis,s_atoms,sup_size)
pw_s.control['prefix'] = pw.control['prefix'][:-1]+"_s'"
if 'celldm(1)' in pw.system: pw_s.system['celldm(1)'] = R[0]*float(pw.system['celldm(1)'])
if 'celldm(2)' in pw.system: pw_s.system['celldm(2)'] = R[1]*float(pw.system['celldm(2)'])/float(pw_s.system['celldm(1)'])
if 'celldm(3)' in pw.system: pw_s.system['celldm(3)'] = R[2]*float(pw.system['celldm(3)'])/float(pw_s.system['celldm(1)'])
if pw.cell_parameters: pw_s.cell_parameters = s_latvec
if 'nbnd' in pw.system: pw_s.system['nbnd'] = sup_size*pw.system['nbnd']
pw_s.system['nat'] = basis*sup_size
pw_s.kpoints = new_kpoints
# Print the output file in the terminal
#print qe
# Write the input file
calculation = ''.join(pw_s.control['calculation'].split('\''))
prefix = ''.join(pw_s.control['prefix'].split('\''))
#pw_s.write('%s.%s'%(prefix,calculation))
pw_s.write('sup_reference.scf')
print 'supercell input file written.'

