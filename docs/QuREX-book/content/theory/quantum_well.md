# Quantum Wells and Walls

## 1. Introduction

We consider a quantum particle of mass {math}` m ` moving in a one-dimensional potential {math}` V(x) `, and study the **stationary states** by solving the **time-independent Schrödinger equation (TISE)**:

<a name='eq:schroedinger'></a>

{math}`
\begin{equation}
\frac{\hbar^2}{2m} \frac{d^2\psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)
\end{equation}
{math}`

The wavefunction {math}` \psi(x) ` must satisfy appropriate **boundary and continuity conditions**, and must be **square-integrable**:

{math}`
\int_{-\infty}^{\infty} |\psi(x)|^2 dx < \infty
{math}`

We study three important cases:

- **Free particle**: {math}` V(x) = 0 `
- **Infinite square well**: {math}` V(x) = \infty ` outside a finite interval
- **Finite square well**: {math}` V(x) = V_0 ` outside a finite interval, with {math}` V_0 < \infty `

---

## 2. Free Particle: {math}` V(x) = 0 `

### Schrödinger Equation:

Set {math}` V(x) = 0 ` in [Schroedinger equation](#eq:schroedinger):
<a name='eq:free_particle_eq'></a>

{math}`
\begin{equation}
\frac{d^2\psi(x)}{dx^2} + k^2 \psi(x) = 0,
\quad \text{with} \quad k = \frac{\sqrt{2mE}}{\hbar}
\end{equation}
{math}`

The general solution is:

{math}`
\psi(x) = Ae^{ikx} + Be^{-ikx}
{math}`

This corresponds to **plane wave solutions**, with energy:

{math}`
\begin{equation}
E = \frac{\hbar^2 k^2}{2m}
\end{equation}
{math}`

These solutions are **not square integrable**, and correspond to **scattering states**. We interpret them as momentum eigenstates:
- {math}` Ae^{ikx} `: right-moving particle with momentum {math}` +\hbar k `
- {math}` Be^{-ikx} `: left-moving particle with momentum {math}` -\hbar k `

---

## 3. Infinite Square Well (Potential Box)

In this case the potential assumes the form:
<a name='eq:inf_well_V'></a>

{math}`
V(x) = \begin{cases}
0 & \text{for } 0 < x < L \\
\infty & \text{otherwise}
\end{cases}
{math}`

- Outside the well ({math}` x \leq 0 ` or {math}` x \geq L `), the wavefunction must vanish:
<a name ='eq:inf_well_BC'></a>

{math}`
\psi(0) = \psi(L) = 0
{math}`

- Inside the well: {math}` V(x) = 0 `

The [TISE](#eq:schroedinger) becomes:

{math}`
\frac{d^2\psi}{dx^2} + k^2 \psi = 0,
\quad k = \frac{\sqrt{2mE}}{\hbar}
{math}`

with general solution:

{math}`
\psi(x) = A\sin(kx) + B\cos(kx)
{math}`

Apply boundary conditions:
- {math}` \psi(0) = 0 \Rightarrow B = 0 `
- {math}` \psi(L) = 0 \Rightarrow \sin(kL) = 0 \Rightarrow k_n = \frac{n\pi}{L} `

The final solutions are then:
<a name='eq:inf_well_solution'></a>

{math}`
\begin{equation}
\psi_n(x) = \sqrt{\frac{2}{L}} \sin\left( \frac{n\pi x}{L} \right), \quad
E_n = \frac{n^2\pi^2\hbar^2}{2mL^2},\quad n = 1,2,3,\dots
\end{equation}
{math}`

#### Key Features:
- Energy spectrum is **discrete** and **non-degenerate**
- No zero-point: {math}` E_1 \neq 0 `
- Energies increase as {math}` n^2 `

---

## 4. Finite Square Well
The potential assumes the form
<a name='eq:finite_well_V'></a>

{math}`
V(x) = \begin{cases}
0 & \text{if } |x| < a \\
V_0 & \text{if } |x| \geq a
\end{cases}
\quad\text{with } V_0 > 0
{math}`

We analyze **bound states** where {math}` E < V_0 `. The [TISE](#eq:schroedinger) must be solved in each region:

---

### Region I: {math}` x < -a `

{math}`
\frac{d^2\psi}{dx^2} = \frac{2m(V_0 - E)}{\hbar^2} \psi = \kappa^2 \psi,
\quad \kappa = \frac{\sqrt{2m(V_0 - E)}}{\hbar}
{math}`

Solution (exponentially decaying):
{math}`
\psi_I(x) = Ae^{\kappa x}, \quad \text{reject } e^{-\kappa x} \text{ (diverges)}
{math}`

---

### Region II: {math}` -a < x < a `

{math}`
\frac{d^2\psi}{dx^2} + k^2 \psi = 0,
\quad k = \frac{\sqrt{2mE}}{\hbar}
{math}`

General solution:

{math}`
\psi_{II}(x) = B\cos(kx) + C\sin(kx)
{math}`

---

### Region III: {math}` x > a `

{math}`
\psi_{III}(x) = De^{-\kappa x}, \quad \text{reject } e^{\kappa x} \text{ (diverges)}
{math}`

---

### Even/Odd Solutions:

We use the symmetry {math}` V(x) = V(-x) \Rightarrow \psi(x) ` is either **even** or **odd**.

#### Even Solutions:

{math}`
\psi(x) = \begin{cases}
A e^{\kappa x} & x < -a \\
B \cos(kx) & |x| < a \\
A e^{-\kappa x} & x > a
\end{cases}
{math}`

Matching at {math}` x = a ` imposes continuity of {math}` \psi ` and {math}` \psi' `:

<a name='eq_even_odd'></a>

{math}`
B\cos(ka) = A e^{-\kappa a},\quad -Bk \sin(ka) = -A \kappa e^{-\kappa a}
\Rightarrow \tan(ka) = \frac{\kappa}{k}
{math}`

#### Odd Solutions:

{math}`
\psi(x) = \begin{cases}
-A e^{\kappa x} & x < -a \\
B \sin(kx) & |x| < a \\
A e^{-\kappa x} & x > a
\end{cases}
{math}`

Matching at {math}` x = a ` gives:

<a name='eq:odd_cond'></a>

{math}`
B \sin(ka) = A e^{-\kappa a},\quad Bk \cos(ka) = -A \kappa e^{-\kappa a}
\Rightarrow -\cot(ka) = \frac{\kappa}{k}
{math}`

---

### Graphical/Numerical Solution:

Equations \eqref{eq:even_cond} and \eqref{eq:odd_cond} are transcendental → solve graphically or numerically for {math}` E `.

Define:

{math}`
z = ka, \quad z_0 = a \sqrt{\frac{2mV_0}{\hbar^2}} \Rightarrow \kappa a = \sqrt{z_0^2 - z^2}
{math}`

Plot left-hand and right-hand sides of:
- {math}` \tan(z) = \sqrt{z_0^2 - z^2} / z ` (even)
- {math}` -\cot(z) = \sqrt{z_0^2 - z^2} / z ` (odd)

---

## 5. Remarks and Limits

- As {math}` V_0 \to \infty `, the equations reduce to the infinite well solution.
- For small {math}` V_0 `, fewer bound states exist.
- The wavefunction **leaks** into classically forbidden regions {math}` |x| > a `: tunneling effect.

---

## 6. Exercises

## 7. References

- D.J. Griffiths, *Introduction to Quantum Mechanics*, Ch. 2–3
- C. Cohen-Tannoudji, *Quantum Mechanics*, Vol. I
- R. Shankar, *Principles of Quantum Mechanics*

