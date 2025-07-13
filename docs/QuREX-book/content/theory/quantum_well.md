# Quantum Wells — "Buca di Potenziale"

## 1. Introduction

We consider a quantum particle of mass \( m \) moving in a one-dimensional potential \( V(x) \), and study the **stationary states** by solving the **time-independent Schrödinger equation (TISE)**:

<a name='eq:schroedinger'></a>

$$
\begin{equation}
\frac{\hbar^2}{2m} \frac{d^2\psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)
\tag{1.1}
\end{equation}
$$

The wavefunction \( \psi(x) \) must satisfy appropriate **boundary and continuity conditions**, and must be **square-integrable**:

$$
\int_{-\infty}^{\infty} |\psi(x)|^2 dx < \infty
$$

We study three important cases:

- **Free particle**: \( V(x) = 0 \)
- **Infinite square well**: \( V(x) = \infty \) outside a finite interval
- **Finite square well**: \( V(x) = V_0 \) outside a finite interval, with \( V_0 < \infty \)

---

## 2. Free Particle: \( V(x) = 0 \)

### Schrödinger Equation:

Set \( V(x) = 0 \) in [Schroedinger equation](#eq:schroedinger):
<a name='eq:free_particle_eq'></a>

$$
\begin{equation}
\frac{d^2\psi(x)}{dx^2} + k^2 \psi(x) = 0,
\quad \text{with} \quad k = \frac{\sqrt{2mE}}{\hbar}
\tag{2.1}
\end{equation}
$$

The general solution is:

$$
\psi(x) = Ae^{ikx} + Be^{-ikx}
$$

This corresponds to **plane wave solutions**, with energy:

$$
\begin{equation}
E = \frac{\hbar^2 k^2}{2m}
\tag{2.2}
\end{equation}
$$

These solutions are **not square integrable**, and correspond to **scattering states**. We interpret them as momentum eigenstates:
- \( Ae^{ikx} \): right-moving particle with momentum \( +\hbar k \)
- \( Be^{-ikx} \): left-moving particle with momentum \( -\hbar k \)

---

## 3. Infinite Square Well (Potential Box)

In this case the potential assumes the form:
<a name='eq:inf_well_V'></a>

$$
V(x) = \begin{cases}
0 & \text{for } 0 < x < L \\
\infty & \text{otherwise}
\end{cases}
\tag{3.1}
$$

- Outside the well (\( x \leq 0 \) or \( x \geq L \)), the wavefunction must vanish:
<a name ='eq:inf_well_BC'></a>

$$
\psi(0) = \psi(L) = 0
\tag{3.2}
$$

- Inside the well: \( V(x) = 0 \)

The [TISE](#eq:schroedinger) becomes:

$$
\frac{d^2\psi}{dx^2} + k^2 \psi = 0,
\quad k = \frac{\sqrt{2mE}}{\hbar}
\tag{3.3}
$$

with general solution:

$$
\psi(x) = A\sin(kx) + B\cos(kx)
$$

Apply boundary conditions:
- \( \psi(0) = 0 \Rightarrow B = 0 \)
- \( \psi(L) = 0 \Rightarrow \sin(kL) = 0 \Rightarrow k_n = \frac{n\pi}{L} \)

The final solutions are then:
<a name='eq:inf_well_solution'></a>

$$
\begin{equation}
\psi_n(x) = \sqrt{\frac{2}{L}} \sin\left( \frac{n\pi x}{L} \right), \quad
E_n = \frac{n^2\pi^2\hbar^2}{2mL^2},\quad n = 1,2,3,\dots
\tag{3.4}
\end{equation}
$$

#### Key Features:
- Energy spectrum is **discrete** and **non-degenerate**
- No zero-point: \( E_1 \neq 0 \)
- Energies increase as \( n^2 \)

---

## 4. Finite Square Well
The potential assumes the form
<a name='eq:finite_well_V'></a>

$$
V(x) = \begin{cases}
0 & \text{if } |x| < a \\
V_0 & \text{if } |x| \geq a
\end{cases}
\quad\text{with } V_0 > 0
\tag{4.1}
$$

We analyze **bound states** where \( E < V_0 \). The [TISE](#eq:schroedinger) must be solved in each region:

---

### Region I: \( x < -a \)

$$
\frac{d^2\psi}{dx^2} = \frac{2m(V_0 - E)}{\hbar^2} \psi = \kappa^2 \psi,
\quad \kappa = \frac{\sqrt{2m(V_0 - E)}}{\hbar}
\tag{4.2}
$$

Solution (exponentially decaying):
$$
\psi_I(x) = Ae^{\kappa x}, \quad \text{reject } e^{-\kappa x} \text{ (diverges)}
$$

---

### Region II: \( -a < x < a \)

$$
\frac{d^2\psi}{dx^2} + k^2 \psi = 0,
\quad k = \frac{\sqrt{2mE}}{\hbar}
\tag{4.3}
$$

General solution:

$$
\psi_{II}(x) = B\cos(kx) + C\sin(kx)
$$

---

### Region III: \( x > a \)

$$
\psi_{III}(x) = De^{-\kappa x}, \quad \text{reject } e^{\kappa x} \text{ (diverges)}
$$

---

### Even/Odd Solutions:

We use the symmetry \( V(x) = V(-x) \Rightarrow \psi(x) \) is either **even** or **odd**.

#### Even Solutions:

$$
\psi(x) = \begin{cases}
A e^{\kappa x} & x < -a \\
B \cos(kx) & |x| < a \\
A e^{-\kappa x} & x > a
\end{cases}
$$

Matching at \( x = a \) imposes continuity of \( \psi \) and \( \psi' \):

<a name='eq_even_odd'></a>

$$
B\cos(ka) = A e^{-\kappa a},\quad -Bk \sin(ka) = -A \kappa e^{-\kappa a}
\Rightarrow \tan(ka) = \frac{\kappa}{k}
\tag{4.4}
$$

#### Odd Solutions:

$$
\psi(x) = \begin{cases}
-A e^{\kappa x} & x < -a \\
B \sin(kx) & |x| < a \\
A e^{-\kappa x} & x > a
\end{cases}
$$

Matching at \( x = a \) gives:

<a name='eq:odd_cond'></a>

$$
B \sin(ka) = A e^{-\kappa a},\quad Bk \cos(ka) = -A \kappa e^{-\kappa a}
\Rightarrow -\cot(ka) = \frac{\kappa}{k}
$$

---

### Graphical/Numerical Solution:

Equations \eqref{eq:even_cond} and \eqref{eq:odd_cond} are transcendental → solve graphically or numerically for \( E \).

Define:

$$
z = ka, \quad z_0 = a \sqrt{\frac{2mV_0}{\hbar^2}} \Rightarrow \kappa a = \sqrt{z_0^2 - z^2}
$$

Plot left-hand and right-hand sides of:
- \( \tan(z) = \sqrt{z_0^2 - z^2} / z \) (even)
- \( -\cot(z) = \sqrt{z_0^2 - z^2} / z \) (odd)

---

## 5. Remarks and Limits

- As \( V_0 \to \infty \), the equations reduce to the infinite well solution.
- For small \( V_0 \), fewer bound states exist.
- The wavefunction **leaks** into classically forbidden regions \( |x| > a \): tunneling effect.

---

## 6. Exercises

## 7. References

- D.J. Griffiths, *Introduction to Quantum Mechanics*, Ch. 2–3
- C. Cohen-Tannoudji, *Quantum Mechanics*, Vol. I
- R. Shankar, *Principles of Quantum Mechanics*

