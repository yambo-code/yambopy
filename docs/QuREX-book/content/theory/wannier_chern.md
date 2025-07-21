# Exciton Chern number 
To compute the exciton Chern number we need to compute overlaps between the periodic part of the exciton wavefunction. In particular, we have to compute the overlap between at point $Q$ belonging to a plane and $Q+\Delta Q_{\hat{i}}/$, where $\Delta Q_{\hat{i}}$ ($ i = 1,2$) denote an increment (of size $1/N_{Q_{i}}$) along the first or second direction in the plane.
This method is referred as the plaquette method.

We want to compute the exciton overlap $M_{\lambda,\lambda^\prime}(\mathbf{Q},\mathbf{Q}+\mathbf{Q_{\hat{i}}})$, in terms of the exciton wavefunction written in the center of mass of the exciton $F_{\lambda,\mathbf{Q}}(\mathbf{R},\mathbf{r})$ {cite}`haber2023maximally` .

$$
\begin{align}
M_{\lambda \lambda^{\prime}}(\mathbf{Q}, \mathbf{Q}+\mathbf{Q_{\hat{i}}})= 
& \langle F_{\lambda,\mathbf{Q}},\mathbf{Q}+\mathbf{Q_{\hat{i}}}\rangle \nonumber \\
& =\sum_{c v \mathbf{k}, c^{\prime} v^{\prime} \mathbf{k}^{\prime}} A_{c v \mathbf{k}+\alpha \mathbf{Q}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}^{\prime}+\alpha \mathbf{Q}+\mathbf{Q_{\hat{i}}}}^{\lambda^{\prime} \mathbf{Q}+\mathbf{Q_{\hat{i}}}} \nonumber\\
&\left[\int_{\mathrm{uc}} d \mathbf{R} \int_{V_{\mathbf{k}}} d \mathbf{r} \chi_{c v \mathbf{k} \mathbf{Q}}^{(\alpha, \beta) \star}(\mathbf{R}, \mathbf{r}) \chi_{c^{\prime} v^{\prime} \mathbf{k}^{\prime} \mathbf{Q}+\mathbf{Q_{\hat{i}}}}^{(\alpha, \mathbf{Q}}(\mathbf{R}, \mathbf{r})\right] \nonumber\\
& =\sum_{c v c^{\prime} v^{\prime} \mathbf{k}} A_{c v \mathbf{k}+\alpha \mathbf{Q}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}^{\prime}+\alpha \mathbf{Q}+\alpha \mathbf{Q_{\hat{i}}}}^{\lambda^{\prime} \mathbf{Q}+\mathbf{Q_{\hat{i}}}}\left\langle u_{c \mathbf{k}+\alpha \mathbf{Q}} \mid u_{c^{\prime} \mathbf{k}+\alpha \mathbf{Q}+\alpha \mathbf{Q_{\hat{i}}}} \right\rangle_{\mathrm{uc}} \nonumber \\
&\left\langle u_{v^{\prime} \mathbf{k}-\beta \mathbf{Q}-\beta \mathbf{Q_{\hat{i}}}} \mid u_{v \mathbf{k}-\beta \mathbf{Q}}\right\rangle_{\mathrm{uc}} \delta_{\mathbf{k} \mathbf{k}^{\prime}} \nonumber\\
& =\sum_{c v c^{\prime} v^{\prime} \mathbf{k}} A_{c v \mathbf{k}}^{\lambda \mathbf{Q} \star} A_{c^{\prime} v^{\prime} \mathbf{k}+\alpha \mathbf{Q_{\hat{i}}}}^{\lambda^{\prime} \mathbf{Q}+\mathbf{Q_{\hat{i}}}}\left\langle u_{c \mathbf{k}} \mid u_{c^{\prime} \mathbf{k}+\alpha \mathbf{Q_{\hat{i}}}}\right\rangle_{\mathrm{uc}}\left\langle u_{v^{\prime} \mathbf{k}-\mathbf{Q}-\beta \mathbf{Q_{\hat{i}}}} \mid u_{v \mathbf{k}-\mathbf{Q}}\right\rangle_{\mathrm{uc}} .
\end{align}
$$ (plaquette-overlap)

<span style="color:darkgreen">
</span>.

# References

```{bibliography}

