---
title: "Improving convergence of Newton's method for degenerate functions"
collection: talks
excerpt: "We explore the application of Newton's method when solving the Butler-Volmer equation, and ways to reduce the number of iterations taken by the solver."
date: 2025-09-20
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

## Introduction
Newton's method is a perpetual workhorse that has proven to be one of the most successful and versatile nonlinear solvers. The framework has been applied to cover a large range of nonlinearity, from s to piecewise-smooth (semismooth) functions [^1,2]. 

### Algorithm 
Let $F \in C^\infty(\mathbb{R}), \; F : \mathbb{R} \rightarrow \mathbb{R}$ be a given smooth function, and we wish to solve $F(x) = 0$. The Newton's method generates a sequence $\{x^{(m)\}$ iteratively: given $x^{(m-1)}$, we obtain $x^(m)$ as

\begin{equation}
\label{eq:Newton_method}
    R^{(m-1)} &=& Fleft(x^{(m-1)} \right),
    \\
    \delta^{(m-1)} &=& - {J^{(m-1)}}^{-1} R^{(m-1)},
    \\
    x^{(m)} &=& x^{(m-1)} + \delta^{(m-1)},
\end{equation}

where $J^{(m-1)} = (F)'(x^{(m-1)}) is the Frechet derivative. 

### Convergence 
The convergence of Newton 

## References
[^1]: C.T. Kelley , *Iterative Methods for Linear and Nonlinear Equations*, 1995, Society for Industrial and Applied Mathematics.
[^2]: Michael Ulbrich, *Semismooth Newton Methods for Variational Inequalities and Constrained Optimization Problems in Function Spaces*, 2011, Mathematical Optimization Society and the Society for Industrial and Applied Mathematics.