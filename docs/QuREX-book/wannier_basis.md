# Wannier representation
We refer to $|u_{n\mathbf{k}}>$ as the Bloch function computed by the DFT code (i.e. QE).
We distinguish between Bloch-like functions in the `Wannier basis`, **Wannier gauge** (which for our purpose can be any localized basis set), and the ones expressed in a Bloch-like band basis **Hamiltonian gauge**.

- **Wannier gauge**: 
  
  $$
  |u^{W}_{n\mathbf{k}}> = \sum_{\mathbf{R}} e^{-i \mathbf{k} \cdot(\hat{\mathbf{r}}-\mathbf{R})}|\mathbf{R} n\rangle
  $$ (eq-wannier-bloch-function)
  
- **Hamiltonian gauge**:
  
  $$
    \left|u_{n \mathbf{k}}^{(\mathrm{H})}\right\rangle=\sum_m\left|u_{m \mathbf{k}}^{(\mathrm{W})}\right\rangle U_{m n}(\mathbf{k})
  $$ (eq-hamiltonian-bloch-function)

where $U$ is the unitary matrix that diagonalizes {eq}`eq-wannier-bloch-function` 

# References

```{bibliography}
