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

where $x \in \mathbb{R}$ is the state variable, $u \in \mathbb{R}$ is the input, $y \in \mathbb{R}$ is the output, $A, B, C, D \in \mathbb{R}$ are constants. The system \ref{eq_state}-\ref{eq_output} may be discretized implicitly or explicitly in time, but that that leads to a time discretization error. Instead, assuming $u$ to be a piecewise constant, i.e., 

\begin{equation}
u(t) = u_k, \; \forall t \in \left(k\tau, (k+1)\tau \right), \; u_k \in \mathbb{R}, \; \forall k \in \mathbb{Z}^{+}, \nonumber
\end{equation}

 a simple integration by parts yields[^1]

\begin{equation}
\label{eq_disc_state}
x_{k+1} = e^{A \tau} x_k + \left(\int_{0}^{k\tau} e^{A z} dz \right) B u_k,
\end{equation}

\begin{equation}
\label{eq_disc_output}
y_k = C x_k + D u_k,
\end{equation}

where $\tau$ is the time step, and $x_k = x(\tau k)$ is the value at the $k^{th}$ time step (in this case the exact value of the continuous solution!), with similar definition for $y_k$.

To understand the discretized system, consider a simple example where $A = D = 0$ and $B = C = 1$. Let $x [m]$ represent the location of an object, and $u [m/s]$ be the velocity which is provided as an input. Then the system \ref{eq_state}-\ref{eq_output} determines the position of the object at any time $t$. In particular, the exact solution is given by $x(t) = \int_{0}^t u(z)dz$, and the output measured is simply the position of the object at $t$, i.e., $y(t) = x(t)$. Using the discretization \ref{eq_disc_state}, we have 

\begin{equation}
\label{eq_disc_motion1}
x_{k+1} = x_k + \tau u_k.
\end{equation}

Equation \ref{eq_disc_motion1} is straightforward and intuitive: it tells us the change in the position of the object in one time step with velocity $u_k$. The corresponding measurement is given using \ref{eq_disc_output}

\begin{equation}
y_k = x_k.
\end{equation}

However, in reality, one might expect some errors to be introduced when measuring the location of the object, or when providing a velocity at a time step. In other words, we expect some noise in the form 

\begin{equation}
\label{eq_disc_motion2}
x_{k+1} = x_k + \tau u_k + w_k
\end{equation}

and 

\begin{equation}
\label{eq_disc_output2}
y_k = x_k + v_k.
\end{equation}

We assume $v_k$ and $w_k$ to be uncorrelated Gaussian random processes with zero mean and given covariance matrices[^2]. The **problem** now becomes as follows: given $\left\{u_k \right\}_k$ and $\left\{y_k \right\}_k$, how can we accurately estimate $\left\{x_k \right\}_k$? 

<br>

We now consider the Kalman filter to estimate the state $\left\{x_k \right\}_k$.

## Kalman filter


## References
[^1]: Chi-Tsong Chen, Linear System Theory and Design, 1999, Oxford University Press, 3rd Edition.
[^2]: Gregory L. Plett, Extended Kalman filtering for battery management systems of LiPB-based HEV battery packs, Part 1, Background, 2004, Journal of Power Sources.
