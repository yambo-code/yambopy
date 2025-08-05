# GW Approximation

## Introduction

The **GW approximation** is a many-body perturbation theory method for calculating quasiparticle energies in materials. It provides accurate band gaps and electronic excitation energies by including exchange-correlation effects beyond density functional theory (DFT).

The name "GW" comes from the Green's function {math}`G` and screened Coulomb interaction {math}`W` that form the basis of the method.

## Theoretical Framework

### Quasiparticle Equation

The quasiparticle energies are solutions to:

```{math}
:label: eq:qp-equation
\left[ \hat{H}_0 + \Sigma(\mathbf{r}, \mathbf{r}', E) \right] \psi_{n\mathbf{k}}(\mathbf{r}) = E_{n\mathbf{k}} \psi_{n\mathbf{k}}(\mathbf{r})
```

where:
- {math}`\hat{H}_0` is the non-interacting Hamiltonian
- {math}`\Sigma` is the self-energy operator
- {math}`\psi_{n\mathbf{k}}` are quasiparticle wavefunctions
- {math}`E_{n\mathbf{k}}` are quasiparticle energies

### Self-Energy in GW

The self-energy in the GW approximation is:

```{math}
:label: eq:gw-self-energy
\Sigma(\mathbf{r}, \mathbf{r}', \omega) = \frac{i}{2\pi} \int d\omega' G(\mathbf{r}, \mathbf{r}', \omega + \omega') W(\mathbf{r}, \mathbf{r}', \omega')
```

### Green's Function

The single-particle Green's function is:

```{math}
:label: eq:greens-function
G(\mathbf{r}, \mathbf{r}', \omega) = \sum_{n\mathbf{k}} \frac{\psi_{n\mathbf{k}}(\mathbf{r}) \psi^*_{n\mathbf{k}}(\mathbf{r}')}{\omega - E_{n\mathbf{k}} + i\eta \text{sgn}(E_F - E_{n\mathbf{k}})}
```

### Screened Coulomb Interaction

The dynamically screened Coulomb interaction is:

```{math}
:label: eq:screened-coulomb
W(\mathbf{r}, \mathbf{r}', \omega) = \int d\mathbf{r}'' \varepsilon^{-1}(\mathbf{r}, \mathbf{r}'', \omega) v(\mathbf{r}'' - \mathbf{r}')
```

where {math}`v(\mathbf{r} - \mathbf{r}')` is the bare Coulomb interaction and {math}`\varepsilon^{-1}` is the inverse dielectric function.

## Practical Implementation

### G₀W₀ Approximation

The most common implementation uses:
- {math}`G_0`: Green's function from DFT eigenvalues and eigenfunctions
- {math}`W_0`: Screened interaction from DFT charge density

The quasiparticle energies are then:

```{math}
:label: eq:g0w0-correction
E_{n\mathbf{k}}^{QP} = E_{n\mathbf{k}}^{DFT} + \langle \psi_{n\mathbf{k}} | \Sigma(E_{n\mathbf{k}}^{QP}) - V_{xc} | \psi_{n\mathbf{k}} \rangle
```

### Computational Steps

1. **DFT Calculation**: Obtain ground state density and wavefunctions
2. **Dielectric Function**: Calculate {math}`\varepsilon(\mathbf{q}, \omega)` 
3. **Screened Interaction**: Compute {math}`W(\mathbf{q}, \omega) = v(\mathbf{q})/\varepsilon(\mathbf{q}, \omega)`
4. **Self-Energy**: Evaluate {math}`\Sigma` matrix elements
5. **Quasiparticle Energies**: Solve for {math}`E_{n\mathbf{k}}^{QP}`

## Connection to BSE

GW provides the quasiparticle energies that serve as input to the [Bethe-Salpeter equation](bse_equation):

- **Single-particle energies**: {math}`\epsilon_{c/v\mathbf{k}} = E_{c/v\mathbf{k}}^{QP}`
- **Band gap correction**: Improved starting point for excitonic calculations
- **Screening consistency**: Same dielectric function used in BSE kernel

## Yambo Implementation

### Key Parameters

- `GbndRnge`: Range of bands for Green's function
- `NGsBlkXp`: G-vectors for polarizability
- `BndsRnXp`: Bands for polarizability calculation
- `EXXRLvcs`: Exchange cutoff
- `PPAPntXp`: Frequency points for plasmon-pole approximation

### Convergence Considerations

Critical parameters for convergence:
- **k-point sampling**: Brillouin zone discretization
- **Band summations**: Number of empty states included
- **G-vector cutoffs**: Plane wave basis size
- **Frequency integration**: Energy mesh for self-energy

## Applications in QuREX

GW results provide:
- Accurate band gaps for [model Hamiltonians](model_hamiltonian)
- Quasiparticle energies for [BSE calculations](bse_equation)
- Screening functions for [Coulomb potentials](coulomb_potential)
- Reference energies for [Wannier interpolation](wannier_basis)

## Limitations and Extensions

### Standard GW Limitations
- Neglects vertex corrections
- Assumes weak correlation
- Single-shot approximation

### Beyond GW
- **Self-consistent GW**: Iterate Green's function and/or screened interaction
- **GW+DMFT**: Include local correlations
- **Vertex corrections**: Include higher-order diagrams
- **GW+BSE**: Combined approach for optical properties

# References

```{bibliography}
```