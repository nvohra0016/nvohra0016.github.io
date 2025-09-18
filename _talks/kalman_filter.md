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

```{math}
:label: my-equation
w_{t+1} = (1 + r_{t+1}) s(w_t) + y_{t+1}
```

See [](#my-equation) for more information!

where $x \in \mathbb{R}$ is the state variable, $y \in \mathbb{R}$ is the output, $A, C \in \mathbb{R}$ are constants. 
The system [(1)](#eq_state)-[(2)](#eq_output) may be discretized implicitly or explicitly in time, but that that leads to a time discretization error. Instead, an integration by 

\begin{eqnarray}
x_{k+1} = A x_k, 
\\
y_k = C x_k.
\end{eqnarray}