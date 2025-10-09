# Bethe-Salpeter Equation (BSE)

## Introduction

The **Bethe-Salpeter Equation (BSE)** provides a rigorous framework within many-body perturbation theory to describe neutral excitations such as excitons, plasmons, and magnons. The BSE incorporates electron-hole interactions beyond the single-particle picture, enabling accurate optical spectra calculations.

The BSE builds upon the **[GW approximation](GW)** for quasiparticle energies, bridging ground-state DFT with excited-state properties across materials from bulk semiconductors to low-dimensional systems.

This chapter follows the ab initio many-body perturbation theory formalism{cite}`sangalli2019many,marini2009yambo` as implemented in Yambo. For implementation details, see the [Yambo BSE flags documentation](../software/yambo/yambo_input_flags.md#bse-bsk).

## Theoretical Foundation

### Many-Body Effects

The transition from DFT to excited-state properties requires addressing:

1. **Quasiparticle Corrections**: Single-particle energy modifications (GW)
2. **Excitonic Effects**: Electron-hole interactions (BSE)

The BSE systematically includes these effects, transforming independent-particle spectra to account for:
- **Bound exciton states** below the gap
- **Continuum resonances** above the gap
- **Oscillator strength redistribution**

## BSE Hamiltonian

In the optical limit ({math}`q = 0`) and Tamm-Dancoff approximation, the BSE Hamiltonian is:

```{math}
:label: eq:bse-hamiltonian
H_{vc\mathbf{k},v'c'\mathbf{k}'} = (\epsilon_{c\mathbf{k}} - \epsilon_{v\mathbf{k}}) \delta_{cc'} \delta_{vv'} \delta_{\mathbf{kk}'} + (f_{c\mathbf{k}} - f_{v\mathbf{k}}) K_{vc\mathbf{k},v'c'\mathbf{k}'}
```

where:
- {math}`\epsilon_{c/v\mathbf{k}}` are quasiparticle energies
- {math}`f_{c/v\mathbf{k}}` are occupation factors
- {math}`K` is the electron-hole interaction kernel

### Interaction Kernel

The kernel decomposes into exchange and correlation parts:

```{math}
:label: eq:kernel
K_{vc\mathbf{k},v'c'\mathbf{k}'} = 2\bar{V}_{vc\mathbf{k},v'c'\mathbf{k}'} - W_{vc\mathbf{k},v'c'\mathbf{k}'}
```

#### Exchange Kernel

The exchange kernel accounts for repulsive electron-hole exchange:

```{math}
:label: eq:exchange-kernel
K^x_{vc\mathbf{k},v'c'\mathbf{k}'} = \bar{V}_{vc\mathbf{k},v'c'\mathbf{k}'} = \frac{1}{\Omega} \sum_{\mathbf{G,G'} \neq \mathbf{0}} v(\mathbf{G}) \langle c\mathbf{k}|e^{i\mathbf{G}\cdot\mathbf{r}}|v\mathbf{k}\rangle \langle v'\mathbf{k}'|e^{-i\mathbf{G}'\cdot\mathbf{r}}|c'\mathbf{k}'\rangle
```

#### Correlation Kernel

The correlation kernel describes screened electron-hole attraction:

```{math}
:label: eq:correlation-kernel
K^c_{vc\mathbf{k},v'c'\mathbf{k}'} = W_{vc\mathbf{k},v'c'\mathbf{k}'} = \frac{1}{\Omega} \sum_{\mathbf{G,G'}} v(\mathbf{q}+\mathbf{G}) \varepsilon^{-1}_{\mathbf{G,G'}}(\mathbf{q}) \langle c\mathbf{k}|e^{i(\mathbf{q}+\mathbf{G})\cdot\mathbf{r}}|c'\mathbf{k}'\rangle \langle v'\mathbf{k}'|e^{-i(\mathbf{q}+\mathbf{G}')\cdot\mathbf{r}}|v\mathbf{k}\rangle
```

## Eigenvalue Problem

The BSE eigenvalue equation is:

```{math}
:label: eq:bse-eigenvalue
\sum_{\mathbf{k}'c'v'} H_{vc\mathbf{k},v'c'\mathbf{k}'} A^{\lambda}_{v'c'\mathbf{k}'} = E_{\lambda} A^{\lambda}_{vc\mathbf{k}}
```

where {math}`A^{\lambda}_{vc\mathbf{k}}` are exciton wavefunctions and {math}`E_{\lambda}` are exciton energies.

## Optical Properties

The macroscopic dielectric function is:

```{math}
:label: eq:dielectric-function
\varepsilon_M(\omega) = 1 - \lim_{\mathbf{q} \to 0} \frac{8\pi}{|\mathbf{q}|^2\Omega} \sum_{vc\mathbf{k}} \sum_{v'c'\mathbf{k}'} \langle v\mathbf{k}|e^{-i\mathbf{q}\cdot\mathbf{r}}|c\mathbf{k}\rangle \langle c'\mathbf{k}'|e^{i\mathbf{q}\cdot\mathbf{r}}|v'\mathbf{k}'\rangle \sum_{\lambda} \frac{A^{\lambda}_{cv\mathbf{k}}(A^{\lambda}_{c'v'\mathbf{k}'})^*}{\omega - E_{\lambda}}
```

## Computational Implementation

### Yambo Parameters

Key Yambo flags for BSE calculations:

- `BSEBands`: Band range for electron-hole pairs
- `BSENGexx`: G-vectors for exchange kernel  
- `BSENGblk`: G-vectors for correlation kernel
- `BandsRnXs`: Bands for screening calculation
- `NGsBlkXs`: G-vectors for screening

### Convergence Considerations

Critical parameters for convergence:
- **k-point sampling**: Brillouin zone discretization
- **Band range**: Number of valence/conduction bands
- **G-vector cutoffs**: Plane wave basis truncation
- **Screening parameters**: Dielectric matrix convergence

## Applications

The BSE enables calculation of:
- Optical absorption spectra
- Exciton binding energies
- Oscillator strengths
- Spatial exciton distributions
- Temperature effects on optical properties

## Connection to QuREX

Within QuREX, BSE results provide:
- Input for [two-particle Hamiltonian](h2p) construction
- Exciton wavefunctions for [Wannier analysis](wannier_exciton)
- Kernel matrix elements for model Hamiltonians
- Benchmark data for [model potentials](coulomb_potential)

# References

```{bibliography}
```