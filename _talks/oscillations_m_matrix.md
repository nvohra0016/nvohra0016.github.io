---
title: "On the Cause of Spurious Oscillations for Stable Numerical Methods"
collection: talks
excerpt: "We explore the cause behind spurious numerical oscillations when solving the heterogeneous heat equation and develop a robust algorithm using appropriate numerical quadrature."
date: 2026-1-23
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

# 1. Introduction

Non-physical oscillations in numerical solutions are almost always linked to the stability of the underlying numerical method. In time dependent problems, the most well known example is the violation of the Courant–Friedrichs–Lewy (CFL) condition for explicit time stepping methods [^1]. The stability for such methods can usually be attained by choosing a small enough time step, which also depends on the spaitial grid size. This however, leads to computationally expensive simulations, requiring more time steps and memory for large time periods. Another way to obtain stability is to use implicit methods, like the backward Euler, which are known to be unconditionally stable. This means that there are no restrictions on the time step size to ensure stability, although a large enough time step comes with loss of accuracy. 

Although it is true that implicit methods do not have a time step restriction, choosing a very small time step can lead to a non-monotone system which also produces oscillations. In this blog post, this is what we precisely demonstrate. We consider the well-known heat equation, and show how oscillations can arise for the stable backward Euler time stepping method. We then propose a remedy for these oscillations, and show its robustness in higher dimensional settings.


# 2. Example: Heat equation

To motivate the existence of spurious oscillations, we consider the heat equation on $\Omega = (0, 1)$

\begin{equation}
\label{eq:heat_eq}
    \partial_t \left(c \theta \right) - \nabla \cdot \left(k \nabla \theta \right) = f, \; \text{ in } \Omega \times (0, T),
\end{equation}

where $\theta$ [$^\circ C$] is the temperature (the primary variable to solve for), $c$ [J/m$^3$ C$^\circ$] is the known volumetric heat capacity of the material, $k = k(x)$ [J/m $^\circ C$ s] is the known thermal conductivity, and $f$ [J/m$^3$ s] is the known external source. Here $T$ [s] denotes the time period. Here we consider homogeneous Dirichlet boundary conditions $\theta(0, t) = \theta(1, t) = 0, \; \forall t \in (0, T)$. 

We discretize our system using piecewise-linear Galerkin elements. In particular, for a given spatial grid size $h > 0$ determined by the number of cells in the grid $N_h$ (i.e., $h = {N_h}^{-1}$), let $V_h$ denote the space of piecewise-linear functions, and let $\\{\phi_{i + \frac{1}{2}}\\}$ denote the basis functions of $V_h$. Following the notation closesly as in our earlier post: [Convergence of Solvers and the Importance of Well-posedness](https://nvohra0016.github.io/talks/elliptic_well_posedness_cg/), we set up the system in a matrix vector form using implicit time stepping as (with $f = 0$)

\begin{equation}
\label{eq:implicit_discretized}
    M \Theta{n} + \tau A \Theta^{n} = M \Theta^{n-1},
\end{equation}

where $\Theta^n \in \mathbb{R}^I$ collects the values of the temperature unknowns at the interior grid points (with $I$ denoting the number of interior grid vertices) at time step $n$, $\tau > 0$ is the time step size, and $M$ and $A$ are the mass and stiffness matrices, respectively. These are defined as

\begin{equation}
    M_{i, j} = \int_0^1 c \phi_i(x) \phi_j(x) dx, \; A_{i, j} = \int_0^1 k(x) \phi_i(x) \phi_j(x) dx.
\end{equation}

It is easy to observe that $M$ and $A$ are symmetric and positive definite (SPD), and have the tri-diagonal $A = \frac{1}{h}$tridiag$(-1, 2, -1)$ and $M = \frac{h}{6}$tridiag$(1,4,1)$ [^2], where $h>0$ is the grid size.

## 2.1. A Stability Estimate of the Implicit Euler Scheme

Let the norm $\\|\cdot \\|_M$ be defined as $\|\U \\|_M = \sqrt{U^T M U}, \; \forall U \in \mathbb{R}^I$. Since $M$ is SPD, it is easy to verify that $\\| \cdot \\|_M$ is a norm. We prove the following stability estimate. 

**Lemma 1.** Let $\Theta^0$ be given. Then, for the scheme given by \ref{equ:implicit_discretized}, we have

\begin{equation}
\label{eq:theorem_stability}
    \\|\Theta^n \\|_M \leq \\|\Theta^{0} \\|_M, \; \forall n \geq 1.
\end{equation}



Note that the estimate \ref{eq:theorem_stability} is independent of $\tau > 0$. We can further use the properties of $M$ to get the estimates...? Let $(U, V) = V^T U, \; \forall U, V \in \mathbb{R}^I$ denote the Euclidean inner product.

\begin{equation}
    \lambda_{min}(M) \\|U \\|_2 \leq \lambda_{max} \|U \|_2.
\end{equation}

 which essentially guarantees that our numerical solution does not blow up as the time step progressess.

*Proof.* Taking the inner product of \ref{eq:implicit_discretized} with $\Theta^{n}$, we have

\begin{equation}
\label{eq:proof1}
    (M\Theta^n - M \Theta^{n-1}, \Theta^n) + \tau (A\Theta^n, \Theta^n) = 0.
\end{equation}

Noting that since $M$ is symmetric, we have $(M \Theta^n, \Theta^{n-1}) = (M \Theta^{n-1}, \Theta^n)$, and hence we have the identity

\begin{equation}
\label{eq:proof2}
    \left(M(\Theta^n - \Theta^{n-1}), \Theta^n \right) = \frac{1}{2} (M\theta^n, \Theta^n) - \frac{1}{2} (M \Theta^{n-1}, \Theta^{n-1}) + \frac{1}{2} (M(\Theta^n - \Theta^{n-1}), \Theta^n - \Theta^{n-1}).
\end{equation}

Substituting \ref{eq:proof2} into \ref{eq:proof1} gives

\begin{equation}
\label{eq:proof3}
    \frac{1}{2} (M\Theta^n, \Theta^n) + \frac{1}{2} (M(\Theta^n - \Theta^{n-1}), \Theta^n - \Theta^{n-1}) + \tau (A\Theta^n, \Theta^n) = \frac{1}{2} (M \Theta^{n-1}, \Theta^{n-1}).
\end{equation}

Now, note that since $M$ and $A$ are SPD, all terms on the LHS of \ref{eq:proof3} are positive, and hence

\begin{equation}
    \frac{1}{2}(M\Theta^n, \Theta^n) \leq \frac{1}{2} (M\Theta^{n-1}, \Theta^{n-1}).
\end{equation}

By induction, we have our result. 

<p style="text-align: right;">&#x25A1;</p>

**Note on the choice of norm.** *If the matrix $M$ were a diagonal matrix, such as $M = h I$ (as we will obtain later), then the norm $\\| \cdot \\|_M = \\| \cdot \\|_2$.*

# 2.2. Test Scenario

We now consider the cooling of a 1D object $1$ [m] long the middle portion of which is kept at $1$ [$^\circ$C], and the rest at $0$ [$^\circ$ C]. That is, we consider the initial condition

$$
\theta_0(x) = 
\begin{cases}
1; & \forall x \in [0.4, 0.6], 
\\
0; & \text{ otherwise}.
\end{cases}
$$

The initial condition is shown in Fig. 1. 

<div align="center">
<img src='/images/m_matrix_oscillations/initial_condition_homogeneous.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 1. Plot showing the initial condition profile considered in the example.
</div>

<br>

We consider material properties similar to water, i.e., we consider $c = 10^6$ [J/ m$^3$ $^\circ$ C] $k = 0.5$ [J/m $^\circ$ C s]. The boundary conditions are as above, i.e., homogeneous Dirichlet. The results for grid size $h = 0.02$ [m] and time step $\tau = 1800$ [hr] are shown in Fig. 2.

<div class='wrapper' align='center'>
<section>
<img id='gif-click' src='/images/neural_network1_blog/training_gif.gif'  width='470' height='470'/>
</section>
</div>

<div align = "center">
 Figure 2. PLACEHOLDER Simulation showing the temperature profile evolution over time.
</div>

<br>

It can be observed that the temperature profile slowly and smoothly decays as the time progresses, and finally reaches a steady state value of $\theta = 0$. 

We now test for the robustness of our algorithm. For large time steps, similar temperature profiles as in Fig. 2. are obtained and nothing interesting happens. On the other hand, let us decrease the time step to gauge the temperature profile over smaller time periods. We consider $\tau = 1$ [s], and re-run the simulation. The results after a few time steps are shown in Fig. 3. 

<div align="center">
<img src='/images/m_matrix_oscillations/oscillations_60_seconds_homogeneous.png' width='380' height='380'>
<img src='/images/m_matrix_oscillations/oscillations_3000_seconds_homogeneous.png' width='380' height='380'>
</div>

<div align = "center">
Figure 3. Results showing the temperature profile near the start of the simulation at $t = 60$ [s] (left) and after consider time steps at $t = 3000$ [s] (right). Notice the oscillations near $x = 0.4$ and $x = 0.6$ when $t = 60$. 
</div>

<br>

It can be observed that now for the first few time steps, spurious oscillations arise in the temperature profile as shown in Fig. 3. (left). The oscillations eventually die out as time progresses, and we get a smooth solution as in the case of $\tau = 3600$. If we look at the $M$-norm values of the temperature, $\\| \theta\\|_M$, over time, we get a monotonically decreasing curve, as expected from Lemma 1; see Fig. 4. 

<div align="center">
<img src='/images/m_matrix_oscillations/norm_homogeneous_small_time_step.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 4. Plot showing the values $\\| \theta \\|_M$ over time.
</div>

<br>

This "overshooting" and "undershooting" behaviour of the function



## References
[^1]: Randall J. LeVeque, *Finite Difference Methods for Ordinary and Partial Differential Equations*, SIAM, 2007.
[^2]: Alexandre Ern, Jean-Luc Guermond, *Theory and Practice of Finite Elements*, 2004, Springer.
[^3]: Fuzhen Zhang, *Matrix Theory: Basic Results and Techniques*, Springer, 2010, Second edition.