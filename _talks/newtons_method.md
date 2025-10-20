---
title: "Accelerating Convergence of Newton's Method for Degenerate Functions"
collection: talks
excerpt: "We explore the application of Newton's method when solving the Butler-Volmer equation, and ways to reduce the number of iterations taken by the solver."
date: 2025-10-19
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

## Introduction to Newton's method
Newton's method is a perpetual workhorse that has proven to be one of the most successful and versatile nonlinear solvers. The framework has been applied to cover a large range of nonlinearity, from s to piecewise-smooth (semismooth) functions [^1,2]. 

### Algorithm 
Let $F \in C^\infty(\mathbb{R}), \; F : \mathbb{R} \rightarrow \mathbb{R}$ be a given smooth function, and we wish to solve $F(x) = 0$. The Newton's method generates a sequence $\\{ x^{(m)} \\}$ iteratively: given $x^{(m-1)}$, we obtain $x^{(m)}$ as

\begin{equation}
\label{eq:Newton_method1}
    R^{(m-1)} = F \left(x^{(m-1)} \right),
\end{equation}

\begin{equation}
\label{eq:Newton_method2}
    \delta^{(m-1)} = - {J^{(m-1)}}^{-1} R^{(m-1)},
\end{equation}

\begin{equation}
\label{eq:Newton_method3}
    x^{(m)} = x^{(m-1)} + \delta^{(m-1)},
\end{equation}

where $J^{(m-1)} = (F)'\left(x^{(m-1)} \right)$ is the Fre´chet derivative (the Jacobian of $F$), and $x^{(0)} = x_0$ is the inital guess which we are given. We iterate till a prescribed tolerance $\epsilon_{tol} > 0$ is reached, i.e., we stop the alogirhm when till $|R^{(m)}| < \epsilon_{tol}$.

**Note on the choice of Jacobian.** In the semismooth framework, the Jacobian $J \in \partial_B F(x)$ is the Clarke's generalized Jacobian which is computed using the B-subdifferential 

\begin{equation}
\partial_B F(x) = \\{ J_F \in \mathbb{R} | \exists \{x_k\} \in D_F, \; x_k \rightarrow x, \; \left(F \right)'(x_k) \rightarrow J_F \\},
\end{equation}

 where for each $x \in D_F \subset \mathbb{R}$, the Fre´chet derivative $F'(x)$ exists [^2].
 

### Convergence of Newton's method
Under the assumptions that $F$ is Lipschitz with bounded derivative, the convergence is well-established for an appropriate initial guess.


## Application to Butler-Volmer equation

The Butler-Volmer equation is used to model the electrochemical reaction kinetics taking place in electrochemical batteries. In particular, they are used to model the lithium ion exchange between the electrode and electrolyte at the microscopic and macroscopic scale[^3]. The equation describes the relationship between the current density $j$ [A/m$^2$] and the electric potential $\phi$ [V] as[^4] 

\begin{equation}
\label{eq:BV}
    j(\phi) = j_0 \left(e^{\frac{\alpha_a z F}{RT} \phi} -  e^{-\frac{\alpha_c z F}{RT} \phi}\right),
\end{equation}

where $j_0$ [A/m$^2$] is the  $\alpha_a$ and $\alpha_c$ [-] are the cathodic and anodic charge transfer coefficients such that $\alpha_a + \alpha_c = 1$, $F \approx 9.648 \times 10^4$ [C/mol] is Faraday's constant, $R \approx 8.314$ [J/K mol] is the gas constant, $z$ is the number of electrons involved in the electrode reaction (for ex. $z = 1$), and $T$ [K] is the temperature. Fig. 1 shows a plot of the current density as a function of the voltage for commonly chosen physical parameters.

<div align="center">
<img src='/images/Newtons_method_images/BV.png' width='380' height='380'>
<img src='/images/Newtons_method_images/BV1.png' width='380' height='380'>
</div>

<div align = "center">
 Figure 1. Plot showing the current density as a function of the potential using the Butler-Volmer equation. Here $\alpha_c = \alpha_a = 0.5$. Left: for $\phi \in [-0.1, 0.1]$. Right: for $\phi \in [-0.5, 0.5].
</div>

<br>

### Example: solving for current
We now consider solving the 1D equation

\begin{equation}
\label{eq:example_bm}
    j(\phi) + \sigma \phi = c,
\end{equation}

where $\sigma$ [S/m] is the conductivity, and the current density $j$ is described by \ref{eq:BV}. For the purposes of this example, we consider $\sigma = 10^6$ to represent the physical values of most materials[^5]. The value of $c$ in the RHS of \ref{eq:example_bm} will be chosen so as to make the solution $\phi$ of \ref{eq:example_bm} closer towards $1$ (which increases the gradient of $j$!). The tolerance is set to $\epsilon_{tol} = 10^{-6}$. 

We test with multiple values $c \in \\{10^4, 10^5, 10^6 \\}$. The absolute value of the residuals $R = j(\phi) +\sigma \phi - c$ is plotted in Fig. 2 below. 

<div align="center">
<img src='/images/Newtons_method_images/alpha_res_smooth.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 2. Plot showing the residuals for different values of $c$. 
</div>

<br>

It can be observed that increasing the value of $c$ increases the number of iterations. This happens because increasing $c$ increases the potential $\phi$, which further increases the gradient $J$, leading to nondegenerate behaviour.  

**Note.** Equations like \ref{eq:example_bm} are encountered as boundary conditions after spatially discretizing the governing equations using, for example, the finite element or finite volume method. Presently, we do not discuss any numerical discretizations, but will revisit that in a future blog post!




## Improving convergence of Newton's method


## References
[^1]: C.T. Kelley, *Iterative Methods for Linear and Nonlinear Equations*, 1995, Society for Industrial and Applied Mathematics.
[^2]: Michael Ulbrich, *Semismooth Newton Methods for Variational Inequalities and Constrained Optimization Problems in Function Spaces*, 2011, Mathematical Optimization Society and the Society for Industrial and Applied Mathematics.
[^3]: F. Brosa et al., *A continuum of physics-based lithium-ion battery models reviewed*, 2022, Progress in Energy, 4.
[^4]: Gregory L. Plett, *Battery Management Systems, Volume 1: Battery Modeling Battery Modeling*, 2015, Artech House Publishers.
[^5]: *Electrical Conductivity of some Common Materials*, Engineering Toolbox, retreived in 2025.