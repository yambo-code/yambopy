# Exciton wavefunction - reference frame and periodic part
The exciton wavefunction $\Psi_{\lambda\mathbf{Q}} (\mathbf{r_e},\mathbf{r_h})$ can be expressed as sum over noninteracting electron-hole products:

$$
\Psi_{\lambda \mathbf{Q}}\left(\mathbf{r}_e, \mathbf{r}_h\right)=\sum_{c v \mathbf{k}} A_{c v \mathbf{k}}^{\lambda \mathbf{Q}} \psi_{c \mathbf{k}}\left(\mathbf{r}_e\right) \psi_{v \mathbf{k}-\mathbf{Q}}^{\star}\left(\mathbf{r}_h\right)
$$ (exc-wavefunction-electronframe)

where $\Psi_{n\mathbf{k}} = e^{i\mathbf{k}\cdot{\mathbf{r}}} u_{n\mathbf{k}}(\mathbf{r})$ denotes a single-particle Bloch state with band index $n$ and crystal momentum $\mathbf{k}$.

We want to compute the exciton overlap $M_{\lambda,\lambda^\prime}(\mathbf{Q},\mathbf{B})$, in terms of the exciton wavefunction written in the center of mass of the exciton $F_{\lambda,\mathbf{Q}}(\mathbf{R},\mathbf{r})$ {cite}`haber2023maximally` .

$$
\begin{align}
M_{\lambda \lambda^{\prime}}(\mathbf{Q}, \mathbf{B})= 
& \langle F_{\lambda,\mathbf{Q}},F_{\lambda,\mathbf{Q+B}}\rangle \nonumber \nonumber\\
& =\sum_{c v \mathbf{k}, c^{\prime} v^{\prime} \mathbf{k}^{\prime}} A_{c v \mathbf{k}+\alpha \mathbf{Q}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}^{\prime}+\alpha \mathbf{Q}+\alpha \mathbf{B}}^{\lambda^{\prime} \mathbf{Q}+\mathbf{B}} \nonumber\\
&\left[\int_{\mathrm{uc}} d \mathbf{R} \int_{V_{\mathbf{k}}} d \mathbf{r} \chi_{c v \mathbf{k} \mathbf{Q}}^{(\alpha, \beta) \star}(\mathbf{R}, \mathbf{r}) \chi_{c^{\prime} v^{\prime} \mathbf{k}^{\prime} \mathbf{Q}+\mathbf{B}}^{(\alpha, \mathbf{Q}}(\mathbf{R}, \mathbf{r})\right] \nonumber\\
& =\sum_{c v c^{\prime} v^{\prime} \mathbf{k}} A_{c v \mathbf{k}+\alpha \mathbf{Q}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}^{\prime}+\alpha \mathbf{Q}+\alpha \mathbf{B}}^{\lambda^{\prime} \mathbf{Q}+\mathbf{B}} \nonumber\\
&\left\langle u_{c \mathbf{k}+\alpha \mathbf{Q}} \mid u_{c^{\prime} \mathbf{k}+\alpha \mathbf{Q}+\alpha \mathbf{B}}\right\rangle_{\mathrm{uc}}\left\langle u_{v^{\prime} \mathbf{k}-\beta \mathbf{Q}-\beta \mathbf{B}} \mid u_{v \mathbf{k}-\beta \mathbf{Q}}\right\rangle_{\mathrm{uc}} \delta_{\mathbf{k} \mathbf{k}^{\prime}} \nonumber\\
& =\sum_{c v c^{\prime} v^{\prime} \mathbf{k}} A_{c v \mathbf{k}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}+\alpha \mathbf{B}}^{\lambda^{\prime},\mathbf{Q}
+\mathbf{B}}\left\langle u_{c \mathbf{k}} \mid u_{c^{\prime} \mathbf{k}+\alpha \mathbf{B}}\right\rangle_{\mathrm{uc}}\left\langle u_{v^{\prime} \mathbf{k}-\mathbf{Q}-\beta \mathbf{B}} \mid u_{v \mathbf{k}-\mathbf{Q}}\right\rangle_{\mathrm{uc}}
\end{align}
$$ (exc-overlap-general)

which for $\alpha = \beta = 1/2$ becomes:

$$
\begin{align}
M_{\lambda \lambda^{\prime}}(\mathbf{Q}, \mathbf{B})= 
& \langle F_{\lambda,\mathbf{Q}},F_{\lambda,\mathbf{Q+B}}\rangle \nonumber \\
= & \sum_{c c^{\prime} v v^{\prime} \mathbf{k}}\left[A_{c v \mathbf{k}}^{\lambda \mathbf{Q}}\right]^{\star} A_{c^{\prime} v^{\prime} \mathbf{k}+\mathbf{B} / 2}^{\lambda^{\prime} \mathbf{Q}+\mathbf{B}} \\
\nonumber
& \times \left \langle u_{c \mathbf{k}} \mid u_{c^{\prime} \mathbf{k}+\mathbf{B} / 2}\right\rangle_{\mathrm{uc}} \left \langle u_{v^{\prime} \mathbf{k}-\mathbf{Q}-\mathbf{B} / 2} \mid u_{v \mathbf{k}-\mathbf{Q}}\right\rangle_{\mathrm{uc}}
\end{align}
$$ (exc-overlap)

The exciton wavefunction in the center of mass reference frame is given by
$$
\begin{align}
F_{\lambda,\mathbf{Q}}
= & \frac{1}{\sqrt{N_{\mathbf{k}}}} \sum_{c v \mathbf{k}} A_{c v \mathbf{k}+\mathbf{Q} / 2}^{\lambda \mathbf{Q}} e^{i \mathbf{k} \cdot \mathbf{r}} \\
& \times u_{c \mathbf{k}+\mathbf{Q} / 2}(\mathbf{R}+\mathbf{r} / 2) u_{v \mathbf{k}-\mathbf{Q} / 2}^{\star}(\mathbf{R}-\mathbf{r} / 2)
\end{align}
$$ (exc-wf-periodic)

# References

```{bibliography}

