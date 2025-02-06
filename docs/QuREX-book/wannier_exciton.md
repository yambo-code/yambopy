# Exciton wavefunction - reference frame and periodic part
The exciton wavefunction $\Psi_{\lambda\mathbf{Q}} (\mathbf{r_e},\mathbf{r_h})$ can be expressed as sum over noninteracting electron-hole products:

$$
\Psi_{S \mathbf{Q}}\left(\mathbf{r}_e, \mathbf{r}_h\right)=\sum_{c v \mathbf{k}} A_{c v \mathbf{k}}^{\lambda \mathbf{Q}} \psi_{c \mathbf{k}}\left(\mathbf{r}_e\right) \psi_{v \mathbf{k}-\mathbf{Q}}^{\star}\left(\mathbf{r}_h\right)
$$ (exc-wavefunction-electronframe)

where $\Psi_{n\mathbf{k}} = e^{i\mathbf{k}\cdot{\mathbf{r}}} u_{n\mathbf{k}}(\mathbf{r})$ denotes a single-particle Bloch state with band index $n$ and crystal momentum $\mathbf{k}$.

We want to compute the exciton overlap $M_{\lambda,\lambda^\prime}(\mathbf{Q},\mathbf{B})$, in terms of the exciton wavefunction written in the center of mass of the exciton $F_{\lambda,\mathbf{Q}}(\mathbf{R},\mathbf{r})$ {cite}`haber2023maximally` .

$$
\begin{align}
M_{S S^{\prime}}(\mathbf{Q}, \mathbf{B})= 
& \langle F_{\lambda,\mathbf{Q}},F_{\lambda,\mathbf{Q+B}}\rangle \nonumber \\
=& \sum_{c c^{\prime} v v^{\prime} \mathbf{k}}\left[A_{c v \mathbf{k}}^{S \mathbf{Q}}\right]^{\star} A_{c^{\prime} v^{\prime} \mathbf{k}+\mathbf{B} / 2}^{S^{\prime} \mathbf{Q}+\mathbf{B}} \\
\nonumber
& \times \left \langle u_{c \mathbf{k}} \mid u_{c^{\prime} \mathbf{k}+\mathbf{B} / 2}\right\rangle_{\mathrm{uc}} \left \langle u_{v^{\prime} \mathbf{k}-\mathbf{Q}-\mathbf{B} / 2} \mid u_{v \mathbf{k}-\mathbf{Q}}\right\rangle_{\mathrm{uc}}
\end{align}
$$ (exc-overlap)

The exciton wavefunction in the center of mass reference frame is given by
$$
\begin{equation}
\begin{aligned}
F_{\lambda,\mathbf{Q}}
= & \frac{1}{\sqrt{N_{\mathbf{k}}}} \sum_{c v \mathbf{k}} A_{c v \mathbf{k}+\mathbf{Q} / 2}^{\lambda \mathbf{Q}} e^{i \mathbf{k} \cdot \mathbf{r}} \\
& \times u_{c \mathbf{k}+\mathbf{Q} / 2}(\mathbf{R}+\mathbf{r} / 2) u_{v \mathbf{k}-\mathbf{Q} / 2}^{\star}(\mathbf{R}-\mathbf{r} / 2)
\end{aligned}
\end{equation}
$$
<span style="color:darkgreen">
How do we make sure that we are in the center of mass reference frame?
From the paper by Haber it looks like we can start from quantities solved in the $vk-Q$ ref. frame
</span>.

# References

```{bibliography}

