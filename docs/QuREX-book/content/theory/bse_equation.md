# Bethe-Salpeter Equation (BSE)

## Introduction to the BSE Formalism

The **Bethe-Salpeter Equation (BSE)** represents the cornerstone of modern theoretical spectroscopy, providing a rigorous framework within **many-body perturbation theory (MBPT)** to describe **neutral excitations** such as excitons, plasmons, and magnons. The BSE extends beyond the single-particle picture by incorporating crucial electron-hole interactions, making it indispensable for calculating **optical spectra** that accurately reproduce experimental observations.

The BSE formalism builds upon the **[GW approximation](GW)** for quasiparticle energies, creating a comprehensive theoretical framework that bridges ground-state DFT calculations with excited-state properties. This approach has proven essential for understanding optical properties in a wide range of materials, from bulk semiconductors to low-dimensional systems.

This chapter follows the **ab initio many-body perturbation theory (AI-MBPT) formalism**{cite}`sangalli2019many,marini2009yambo` as implemented in the **Yambo code**. For practical implementation details and input parameters, refer to the [Yambo BSE flags documentation](../software/yambo/yambo_input_flags.md#bse-bsk).

## Theoretical Foundation

### From Independent Particles to Correlated Excitations

The transition from ground-state DFT to excited-state properties requires addressing two fundamental many-body effects:

1. **Quasiparticle Corrections**: Modification of single-particle energies due to exchange-correlation effects (addressed by GW)
2. **Excitonic Effects**: Electron-hole interactions that can dramatically modify optical spectra (addressed by BSE)

The BSE provides a systematic way to include these correlated effects, transforming the independent-particle absorption spectrum into a realistic description that accounts for:
- **Bound exciton states** below the quasiparticle gap
- **Continuum resonances** above the gap
- **Oscillator strength redistribution** due to many-body effects



### Connection to Experimental Observables
- In the **optical limit** \( q = 0 \), we compute the system’s **response to light absorption**.
- Solving the BSE provides **exciton energies and wavefunctions**, enabling **optical spectra** calculations, which can be directly compared to experimental **absorption spectra** and **electron energy loss spectra (EELS)**.

In the **Tamm-Dancoff approximation (TDA)**, including **local field effects** {cite}`onida2002electronic`, the **two-particle BSE Hamiltonian** for a **non-spin-polarized system** in the **optical limit** \( q = 0 \) is given by:
`yambo -o b -k sex -b`

```{math}
H_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}=
\left(\epsilon_{c \mathbf{k}}-\epsilon_{v \mathbf{k}}\right) \delta_{c, c^{\prime}} \delta_{v, v^{\prime}} \delta_{\mathbf{k k}^{\prime}} 
+ \left(f_{c \mathbf{k}}-f_{v \mathbf{k}}\right)\left[2 \bar{V}_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}} - 
 W_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}\right]
```

### Physical Meaning of Terms
- The indices **vc𝑘** denote the **valence and conduction bands** at **quasiparticle momentum** $ \mathbf{k} $. \
    Yambo: `% BSEBands lower band | upper band |` Bands range.
- The term $ \left(\epsilon_{c \mathbf{k}}-\epsilon_{v \mathbf{k}}\right) $ represents the **quasiparticle energy differences**.\
    Yambo: `KfnQPdb=" E < ./SAVE/ndb.QP"` Location of QP corrections database from previous GW calculation.\
    `% KfnQP_E scissor | stretch conduction | stretch valence ` QP corrections parameter. 
- The second term, 
  $\left[2 \bar{V}_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}} - W_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}\right]$
  forms the **kernel** $ K $, where the key **electron-hole interactions** occur.

For clarity, we decompose the **kernel** into:
1. The **exchange part**: $ K^x $ (Hartree potential contribution).
2. The **correlation part**: $ K^c $ (electron-hole attraction).

### Exchange and Correlation Kernels
The **exchange kernel** $ K^x $ accounts for the **repulsive electron-hole exchange interaction**, given by:

```{math}
K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^x=
\bar{V}_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}=
\frac{1}{\Omega} \sum_{\mathbf{G, G^{\prime}} \neq \mathbf{0}} v(\mathbf{G})
\left\langle c \mathbf{k}\left|e^{i \mathbf{G r}}\right| v \mathbf{k}\right\rangle
\left\langle v^{\prime} \mathbf{k}^{\prime}\left|e^{-i \mathbf{G}^{\prime} \mathbf{r}}\right| c^{\prime} \mathbf{k}^{\prime}\right\rangle
```

 - Yambo: `BSENGexx = ` G Compontents of Hartree potential included in the summation. 

The **correlation kernel** $ K^c $ describes the **screened electron-hole attraction**, incorporating many-body effects via the **inverse dielectric function** $ \varepsilon^{-1}$:

```{math}
K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^c=
W_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}=
\frac{1}{\Omega} \sum_{\mathbf{G}, \mathbf{G}^{\prime}} v(\mathbf{q}+\mathbf{G}) 
\varepsilon_{\mathbf{G}, \mathbf{G}^{\prime}}^{-1}(\mathbf{q})
\left\langle c \mathbf{k}\left|e^{i(\mathbf{q}+\mathbf{G}) \mathbf{r}}\right| c^{\prime} \mathbf{k}^{\prime}\right\rangle 
\left\langle v^{\prime} \mathbf{k}^{\prime}\left|e^{-i\left(\mathbf{q}+\mathbf{G}^{\prime}\right) \mathbf{r}}\right| v \mathbf{k}\right\rangle 
\delta_{\mathbf{q k}-\mathbf{k}^{\prime}}
```
- Yambo: `BSENGblk = `G included in sum for the screened interaction 


The **inverse dielectric screening matrix** $ \varepsilon_{\mathbf{G}, \mathbf{G}^{\prime}}^{-1}(\mathbf{q})$ in $ K^c $ is computed using the **Yambo flag `em1s`**.
This dielectric screening is also a crucial component of the **GW approximation**, linking the BSE formalism to quasiparticle corrections.
- Yambo: `%BandsRnXs lower bound | upper bound |` bands included in screening. \
 `NGsBlkXs = ` Gs included in screening. \
 `LongDrXS x | y | z | ` Compontents of E-field 
   



## Eigenvalue problem

The BSE can be recast into an eigenvalue equation 

```{math}
:label:eq-BSE
\left(\varepsilon_{c\mathbf{k}}
-\varepsilon_{\mathrm{v}\mathbf{k-Q}}\right)
A_{\mathrm{vc}\mathbf{k}}^\lambda+\sum_{\mathbf{k}^{\prime}
c^{\prime} \mathbf{v}^{\prime}}
K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{\mathbf{Q}}
A_{v^{\prime} c^{\prime} \mathbf{k}^{\prime}}^{\lambda,\mathbf{Q}}
=E_{\lambda,\mathbf{Q}} A_{v c \mathbf{k}}^{\lambda,\mathbf{Q}}
```

 $A^{\lambda}_{v c\mathbf{k}}$ are the eigenvectors, $E_{\lambda,\mathbf{Q}}$ are 
the energy of excitonic transition $\lambda$, and $\mathbf{Q}$ the momentum transfer between an electron at $\mathbf{k}$ and a hole at $\mathbf{k-Q}$.
 - Yambo: `% BLongDir x | y | z | %` Direction of longitudinal perturbation.

## Diagonalisation

The macroscopic dielectric function is obtained with the excitonic eigenvectors $A^{\lambda}$ and eigenenergies $E_{\lambda}$: 

```{math}
\varepsilon_{M}(\omega) = 1 - \text{lim}_{\mathbf{q} \rightarrow 0} \dfrac{8\pi}{|\mathbf{q}|^2\Omega} \sum_{v,c,\mathbf{k}}\sum_{v^{\prime},c^{\prime},\mathbf{k}^{\prime}}
\left\langle v \mathbf{k} - \mathbf{q}\left|e^{-i\mathbf{q r}}\right| c \mathbf{k}\right\rangle
\left\langle c^{\prime} \mathbf{k}^{\prime} - \mathbf{q}\left|e^{i\mathbf{q r}}\right| v^{\prime} \mathbf{k^{\prime}- q}\right\rangle
\sum_{\lambda}\dfrac{A^{\lambda}_{cv\mathbf{k}}(A^{\lambda}_{cv\mathbf{k}})^*}{\omega - E_{\lambda}}
```



The kernel $K$ contains the electron-hole Coulomb interaction matrix elements, which is the sum of the direct $K^{d}$ and exchange $K^{x}$ terms and can be written as

$$
K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{\mathbf{Q}} =     K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{d,\mathbf{Q}} +     K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{x,\mathbf{Q}}
$$

Following the prescription of {cite}`dias2023wantibexos` we can compute the direct and exchange terms as follows:

$$
K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{d,\mathbf{Q}} = V(\mathbf{k}-\mathbf{k}^\prime)\langle c,\mathbf{k}|c^\prime,\mathbf{k}^\prime\rangle\langle v^\prime, \mathbf{k}^\prime-\mathbf{Q}|v\mathbf{k-Q}\rangle \\
K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{x,\mathbf{Q}} = -V(\mathbf{k}-\mathbf{k}^\prime)\langle c,\mathbf{k}| v,\mathbf{k}-\mathbf{Q}\rangle\langle v^\prime, \mathbf{k}^\prime-\mathbf{Q}|c^\prime\mathbf{k}^\prime-\mathbf{Q}\rangle \\    
$$ (eq:kernels)


