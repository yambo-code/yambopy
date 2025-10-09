# Exciton Chern number 
To compute the exciton Chern number we need to compute overlaps between the periodic part of the exciton wavefunction. In particular, we have to compute the overlap between at point {math}`Q` belonging to a plane and {math}`Q+\Delta Q_{\hat{i}}/`, where {math}`\Delta Q_{\hat{i}}` ({math}` i = 1,2`) denote an increment (of size {math}`1/N_{Q_{i}}`) along the first or second direction in the plane.
This method is referred as the plaquette method.

We want to compute the exciton overlap {math}`M_{\lambda,\lambda^\prime}(\mathbf{Q},\mathbf{Q}+\mathbf{Q_{\hat{i}}})`, in terms of the exciton wavefunction written in the center of mass of the exciton {math}`F_{\lambda,\mathbf{Q}}(\mathbf{R},\mathbf{r})` {cite}`haber2023maximally` .

{math}`
\begin{align}
M_{\lambda \lambda^{\prime}}(\mathbf{Q}, \mathbf{Q}+\mathbf{Q_{\hat{i}}})= 
& \langle F_{\lambda,\mathbf{Q}},\mathbf{Q}+\mathbf{Q_{\hat{i}}}\rangle \nonumber \\
& =\sum_{c v \mathbf{k}, c^{\prime} v^{\prime} \mathbf{k}^{\prime}} A_{c v \mathbf{k}+\alpha \mathbf{Q}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}^{\prime}+\alpha \mathbf{Q}+\mathbf{Q_{\hat{i}}}}^{\lambda^{\prime} \mathbf{Q}+\mathbf{Q_{\hat{i}}}} \nonumber\\
&\left[\int_{\mathrm{uc}} d \mathbf{R} \int_{V_{\mathbf{k}}} d \mathbf{r} \chi_{c v \mathbf{k} \mathbf{Q}}^{(\alpha, \beta) \star}(\mathbf{R}, \mathbf{r}) \chi_{c^{\prime} v^{\prime} \mathbf{k}^{\prime} \mathbf{Q}+\mathbf{Q_{\hat{i}}}}^{(\alpha, \mathbf{Q}}(\mathbf{R}, \mathbf{r})\right] \nonumber\\
& =\sum_{c v c^{\prime} v^{\prime} \mathbf{k}} A_{c v \mathbf{k}+\alpha \mathbf{Q}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}^{\prime}+\alpha \mathbf{Q}+\alpha \mathbf{Q_{\hat{i}}}}^{\lambda^{\prime} \mathbf{Q}+\mathbf{Q_{\hat{i}}}}\left\langle u_{c \mathbf{k}+\alpha \mathbf{Q}} \mid u_{c^{\prime} \mathbf{k}+\alpha \mathbf{Q}+\alpha \mathbf{Q_{\hat{i}}}} \right\rangle_{\mathrm{uc}} \nonumber \\
&\left\langle u_{v^{\prime} \mathbf{k}-\beta \mathbf{Q}-\beta \mathbf{Q_{\hat{i}}}} \mid u_{v \mathbf{k}-\beta \mathbf{Q}}\right\rangle_{\mathrm{uc}} \delta_{\mathbf{k} \mathbf{k}^{\prime}} \nonumber\\
& =\sum_{c v c^{\prime} v^{\prime} \mathbf{k}} A_{c v \mathbf{k}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}+\alpha \mathbf{Q_{\hat{i}}}}^{\lambda^{\prime} \mathbf{Q}+\mathbf{Q_{\hat{i}}}}\left\langle u_{c \mathbf{k}} \mid u_{c^{\prime} \mathbf{k}+\alpha \mathbf{Q_{\hat{i}}}}\right\rangle_{\mathrm{uc}}\left\langle u_{v^{\prime} \mathbf{k}-\mathbf{Q}-\beta \mathbf{Q_{\hat{i}}}} \mid u_{v \mathbf{k}-\mathbf{Q}}\right\rangle_{\mathrm{uc}} .
\end{align}
{math}` (plaquette-overlap)
<span style="color:darkgreen">
</span>.

General formula for computation of Chern number
{math}`
\begin{aligned}
&C_n=\frac{1}{2 \pi} \int_{\mathrm{BZ}} \Omega_n(\mathbf{k}) d^2 k\\
&\Omega_n(\mathbf{k})=i \sum_{m \neq n} \frac{\left\langle u_{n \mathbf{k}}\right| \partial_{k_x} H(\mathbf{k})\left|u_{m \mathbf{k}}\right\rangle\left\langle u_{m \mathbf{k}}\right| \partial_{k_y} H(\mathbf{k})\left|u_{n \mathbf{k}}\right\rangle-(x \leftrightarrow y)}{\left(E_n(\mathbf{k})-E_m(\mathbf{k})\right)^2}
\end{aligned}
{math}` (chern-general)

# References

```{bibliography}

