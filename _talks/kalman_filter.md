---
title: "An Introduction to Kalman Filters"
collection: talks
excerpt: "We present an investigation into Kalman filters and their applications in BMS"
date: 09-2025
---

## Introduction

Consider a simple 1D linear state-space system of the form 

$$
\begin{aligned}
\label{eq_state}
\dot{x(t)} = A x(t),
\\
\label{eq_output}
y(t) = C x(t),
\end{aligned}
$$

where $x \in \mathbb{R}$ is the state variable, $y \in \mathbb{R}$ is the output, $A, C \in \mathbb{R}$ are constants. Discretizing [%s](#eq_state)-[%s](#eq_output) explicitly in time gives

\begin{eqnarray}
x_{k+1} = A x_k, 
\\
y_k = C x_k.
\end{eqnarray}