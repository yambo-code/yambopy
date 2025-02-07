# BSE Equation

The Bethe-Salper Equations (BSE) allows computing neutral excitations for electron-hole interactions.
Following the Ai-MBPT formalism {cite}`sangalli2019many, marini2009yambo` the exciton energies and wavefunctions, are obtained solving the BSE equation in the Tamm-Dancoff approximation including local field effects {cite}`onida2002electronic`. 

The BSE can be recast into an eigenvalue equation 

$$
    \left(\varepsilon_{c\mathbf{k}}^{\mathrm{GW}}
    -\varepsilon_{\mathrm{v}\mathbf{k-Q}}^{\mathrm{GW}}\right)
    A_{\mathrm{vc}\mathbf{k}}^\lambda+\sum_{\mathbf{k}^{\prime}
    c^{\prime} \mathbf{v}^{\prime}}
    K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{\mathbf{Q}}
    A_{v^{\prime} c^{\prime} \mathbf{k}^{\prime}}^{\lambda,\mathbf{Q}}
    =E_{\lambda,\mathbf{Q}} A_{v c \mathbf{k}}^{\lambda,\mathbf{Q}}
$$ (eq-BSE)

where $\varepsilon_{c\mathbf{k}/v\mathbf{k}}$ are quasi-particle band energies, 
 $A^{\lambda}_{v c\mathbf{k}}$ are the BSE coefficients, $E_{\lambda,\mathbf{Q}}$ %are 
the energy of exciton $\lambda$, and $\mathbf{Q}$ the momentum transfer between an electron at $\mathbf{k}$ and a hole at $\mathbf{k-Q}$.
The kernel $K$ contains the electron-hole Coulomb interaction matrix elements, which is the sum of the direct $K^{d}$ and exchange $K^{x}$ terms and can be written as
$$
    K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{\mathbf{Q}} =     K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{d,\mathbf{Q}} +     K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{x,\mathbf{Q}}
$$

Following the prescription of \textit{Dias et al.} {cite}`dias2023wantibexos`} we can compute the direct and exchange terms as follows:

$$
    K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{d,\mathbf{Q}} = V(\mathbf{k}-\mathbf{k}^\prime)\langle
    c,\mathbf{k}|
    c^\prime,\mathbf{k}^\prime\rangle\langle
    v^\prime, \mathbf{k}^\prime-\mathbf{Q}|v\mathbf{k-Q}\rangle \\
    K_{\substack{\mathrm{vc\mathbf{k}} \\ v^{\prime} c^{\prime} \mathbf{k}^{\prime}}}^{x,\mathbf{Q}} = -V(\mathbf{k}-\mathbf{k}^\prime)\langle
    c,\mathbf{k}|
    v,\mathbf{k}-\mathbf{Q}\rangle\langle
    v^\prime, \mathbf{k}^\prime-\mathbf{Q}|c^\prime\mathbf{k}^\prime-\mathbf{Q}\rangle \\    
$$ (eq:kernels)