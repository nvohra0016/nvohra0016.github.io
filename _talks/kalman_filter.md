---
title: "An Introduction to Kalman Filters"
collection: talks
excerpt: "We present an investigation into Kalman filters and their applications in BMS"
date: 09-2025
---

## Introduction

Consider a simple 1D linear state-space system of the form 

\begin{equation}
\label{eq_state}
\dot{x(t)} = A x(t) + B u(t),
\end{equation}
\begin{equation}
\label{eq_output}
y(t) = C x(t) + D u (t), \; \forall t \in (0, T),
\end{equation}

where $x \in \mathbb{R}$ is the state variable, $u \in \mathbb{R}$ is the input, $y \in \mathbb{R}$ is the output, $A, B, C, D \in \mathbb{R}$ are constants. The system \ref{eq_state} \eqref{eq_state} [(1)](#eq_state)-[(2)](#eq_output) may be discretized implicitly or explicitly in time, but that that leads to a time discretization error. Instead, assuming $u$ to be a piecewise constant, a simple integration by parts yields[^1]

\begin{equation}
x_{k+1} = e^{A \tau} x_k + \left(\int_{0}^{k\tau} e^{A z} dz \right) B u_k,
\end{equation}

\begin{equation}
y_k = C x_k + D u_k,
\end{equation}

where $\tau$ is the time step, and $x_k = x(\tau k)$ is the value at the $k^{th}$ time step (in this case the exact value of the continuous solution!). 

To understand the discretized system, consider a simple example where $A = D = 0$ and $B = C = 1$. Let $x [m]$ represent the location of an object, and $u [m/s]$ be the velocity which is provided as an input. Then the system [(1)](#eq_state)-[(2)](#eq_output) determines the position of the object at any time $t$. In particular, the exact solution is given by $x(t) = \int_{0}^t u(z)dz$, and the output measured is simply the position of the object at $t$, i.e., $y(t) = x(t)$. Using the discretization, we have 



## References
[^1]: Chi-Tsong Chen, Linear System Theory and Design, 1999, 3rd Edition.
