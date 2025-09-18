---
title: "An Introduction to Kalman Filters"
collection: talks
excerpt: "We present an investigation into Kalman filters and their applications in BMS"
date: 09-2025
---

## Introduction

Consider a simple 1D linear state-space system of the form 

\begin{eqnarray}
\label{eq_state}
\dot{x(t)} = A x(t) + B u(t),
\end{eqnarray}
\begin{eqnarray}
\label{eq_output}
y(t) = C x(t),
\end{eqnarray}

where $x \in \mathbb{R}$ is the state variable, $u \in \mathbb{R}$ is the input, $y \in \mathbb{R}$ is the output, $A, C \in \mathbb{R}$ are constants. 
The system [(1)](#eq_state)-[(2)](#eq_output) may be discretized implicitly or explicitly in time, but that that leads to a time discretization error. Instead, assuming $u$ to be a piecewise constant, a simple integration by parts yields^[1]

\begin{eqnarray}
x_{k+1} = e^{A \tau} x_k + \left(\int_{0}^{k\tau} e^{A y} dy \right) B u_k,
\end{eqnarray}

\begin{eqnarray}
y_k = C x_k.
\end{eqnarray}

## References
[^1]. Chi-Tsong Chen, Linear System Theory and Design, 1999, 3rd Edition.
