# BSE equation

The BSE Hamiltonian is construcuted in `wann_H2P.py` using the energies and eigenvectors computed in `wann_model`.
The BSE Hamiltonian can be built from:
1) A model {doc}`coulomb_potential`
2) From Yambo Kernel `YamboBSEKernelDB`
3) It can be reconstructured from `YamboExcitonDB` and the already computed by eigenvalues and eigenvectors by Yambo


$$
H^{2P}_{\lambda,\lambda^\prime} =
$$ (eq:H2P)
