---
title: "Improving convergence of Newton's method for degenerate functions"
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
Let $F \in C^\infty(\mathbb{R}), \; F : \mathbb{R} \rightarrow \mathbb{R}$ be a given smooth function, and we wish to solve $F(x) = 0$. The Newton's method generates a sequence $\{x^{(m)\}$ iteratively: given $x^{(m-1)}$, we obtain $x^(m)$ as

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

where $J^{(m-1)} = (F)'\left(x^{(m-1)} \right)$ is the Fre´chet derivative (the Jacobian of $F$), and $x^{(0)} = x_0$ is the inital guess which we are given. 

<br>

**Note** In the semismooth framework, the Jacobian $J \in \partial_B F(x)$ is the Clarke's generalized Jacobian which is computed using the B-subdifferential 

\begin{equation}
\partial_B F(x) = \{ J_F \in \mathbb{R} | \exists \{x_k\} \in D_F, \; x_k \rightarrow x, \; \left(F \right)'(x_k) \rightarrow J_F\},
\end{equation}

 where for each $x \in D_F \subset \mathbb{R}$, the Fre´chet derivative $F'(x)$ exists [^2].
 

### Convergence 
Under the assumptions that $F$ is Lipschitz with bounded derivative, the convergence is well-established for an appropriate initial guess.


## Application to Butler-Volmer equation

The Butler-Volmer equation is well-known in the areas of electrochemical batteries. In particular, the equation computes describes the relationship between the current density $j$ [A/m$^2$] and the voltage $v$ [V] as[^3] 

\begin{equation}
    j(v) = j_0 \left(e^{\frac{\alpha_a z F}{RT} v} -  e^{-\frac{\alpha_c z F}{RT} v}\right),
\end{equation}

where $j_0$ [A/m$^2$] is the  $alpha_a$ and $\alpha_c$ [-] are the cathodic and anodic charge transfer coefficients such that $\alpha_a + \alpha_c = 1$, $F \approx 9.648 \times 10^4$ [C/mol] is Faraday's constant, $R \approx 8.314$ [J/K mol] is the gas constant, $z$ is the number of electrons involved in the electrode reaction (for ex. $z = 1$), and $T$ [K] is the temperature. Fig. 1 shows a plot of the current density as a function of the voltage for commonly chosen physical parameters.

<div align="center">
<img src='/images/Kalman_filter/heat_ex_temperature_results.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 3. Plot showing the estimated temperature (left) and predicted and corrected variance (right).
</div>


## Improving convergence of Newton's method


## References
[^1]: C.T. Kelley, *Iterative Methods for Linear and Nonlinear Equations*, 1995, Society for Industrial and Applied Mathematics.
[^2]: Michael Ulbrich, *Semismooth Newton Methods for Variational Inequalities and Constrained Optimization Problems in Function Spaces*, 2011, Mathematical Optimization Society and the Society for Industrial and Applied Mathematics.
[^3]: Gregory L. Plett, *Battery Management Systems, Volume 1: Battery Modeling Battery Modeling*, 2015, Artech House Publishers.