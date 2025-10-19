---
title: "Improving convergence of Newton's method for degenerate functions"
collection: talks
excerpt: "We explore the application of Newton's method when solving the Butler-Volmer equation, and ways to reduce the number of iterations taken by the solver."
date: 2025-10-19
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

## Introduction
Newton's method is a perpetual workhorse that has proven to be one of the most successful and versatile nonlinear solvers. The framework has been applied to cover a large range of nonlinearity, from s to piecewise-smooth (semismooth) functions [^1,2]. 

### Algorithm 
Let $F \in C^\infty(\mathbb{R}), \; F : \mathbb{R} \rightarrow \mathbb{R}$ be a given smooth function, and we wish to solve $F(x) = 0$. The Newton's method generates a sequence $\{x^{(m)\}$ iteratively: given $x^{(m-1)}$, we obtain $x^(m)$ as

\begin{equation}
\label{eq:Newton_method1}
    R^{(m-1)} = Fleft(x^{(m-1)} \right),
\end{equation}

\begin{equation}
\label{eq:Newton_method2}
    \delta^{(m-1)} = - {J^{(m-1)}}^{-1} R^{(m-1)},
\end{equation}

\begin{equation}
\label{eq:Newton_method3}
    x^{(m)} = x^{(m-1)} + \delta^{(m-1)},
\end{equation}

where $J^{(m-1)} = (F)'\left(x^{(m-1)} \right)$ is the Frechet derivative (the Jacobian of $F$).

**Note** In the semismooth framework, the Jacobian $J \in \partial_B F(x)$ is the Clarke's generalized Jacobian which is computed using the B-subdifferential $\partial_B F(x) = \{J_F \in \mathbb{R} \; | \; \exists \{x_k\} \in D_F, \; x_k \rightarrow x. \; \left(F \right)'(x_k) \rightarrow J_F\}
 

### Convergence 
The convergence of Newton 

## References
[^1]: C.T. Kelley , *Iterative Methods for Linear and Nonlinear Equations*, 1995, Society for Industrial and Applied Mathematics.
[^2]: Michael Ulbrich, *Semismooth Newton Methods for Variational Inequalities and Constrained Optimization Problems in Function Spaces*, 2011, Mathematical Optimization Society and the Society for Industrial and Applied Mathematics.