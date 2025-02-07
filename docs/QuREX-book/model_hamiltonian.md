(model_ham_intro)=
# Introduction 
The workhorse method for tight-binding (TB) parametrization of the electronic properties of a system is the orthogonal MLWF-TB parametrization {cite}`marzari2012maximally`. In this approach, the electronic TB Hamiltonian ($H$) is extracted from a DFT calculation employing, for example, open-source codes such as Quantum Espresso {cite}`giannozzi2009quantum,giannozzi2017advanced` and Wannier90 {cite}`mostofi2008wannier90`.
In this way, one obtains a real space representation of the electronic Hamiltonian $H_{nm}(\mathbf{R})$, where $\mathbf{R}$ are the lattice vectors lying in a supercell conjugate to the $\mathbf{k}$-mesh and $n$,$m$ label the electronic band indices.
Then, one can use via a Slater-Koster interpolation scheme {cite}`yates2007spectral` the reciprocal space Hamiltonian ($H_{nm}(\mathbf{k})$) on a finer mesh of $\mathbf{k}$-points with respect to the one used in the DFT calculation as

$$
H_{n m}\left(\mathbf{k}\right)=\sum_{\mathbf{R}} \mathrm{e}^{\mathrm{i} \mathbf{k} \cdot \mathbf{R}} H_{n m}(\mathbf{R}) 
$$ (HR)

After diagonalization the resulting eigenvectors are in the **Hamilton gauge** {doc}`model_hamiltonian`