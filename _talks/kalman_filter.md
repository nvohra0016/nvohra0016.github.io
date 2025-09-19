---
title: "An Introduction to Kalman Filters"
collection: talks
excerpt: "We present an investigation into Kalman filters and their applications in BMS"
date: 09-2025
---

## Introduction

We start by considering a simple example of heat conduction and temperature change in an object over a time period $(0, T)$. Consider an object with volumetric heat capacity $c$ [J/m$^3$ $^\circ C$]. We add heat to this object at the rate of $f(t)$ [J/m$^3$ s]. Ignoring any spatial variation, the temperature $\theta(t)$ [$^\circ$ C] of the object can be determined using[^1]

\begin{equation}
\label{eq:heat_eq}
 \dot{(c \theta)} = f \text{ in } (0, T).
\end{equation}

Now, let us measure the temperature using an external sensor, and we denote the measurement using $y(t)$ [$^\circ$ C]; see figure below for a schematic. 

<div align="center">
<img src='/images/heat_kalman_filter_example.png' width='250' height='250'>
</div>

Then, the measurement can be represented by

\begin{equation}
\label{eq:measurement_heat}
y(t) = \theta(t), \; \forall t \in (0, T).
\end{equation}

<br>

Let us proceed with a discretization of \ref{eq:heat_eq}. Let $\tau$ denote the given time step and let $\tilde{f} \approx f, \; \tilde{f} \in C^0(0, T)$ be a piecewise-constant approximation given by

\begin{equation}
\label{eq:f_assum}
\tilde{f}(t) = f_k = f(k \tau), \; \forall t \in \left(k\tau, (k+1)\tau \right), \; \forall k \in \mathbb{Z}^{+}.
\end{equation}

Under the assumption \ref{eq:f_assum}, the solution to 

\begin{equation}
\label{eq:heat_eq1}
\dot{c \theta} = \tilde{f} \text{ in } (0, T),
\end{equation}

is given by $\theta(t) = \int_0^t c^{-1} \tilde{f}(z)dz \in H^1(0, T)$, which is piecewise-linear. To estimate $\theta$, we can rewrite \ref{eq:heat_eq1} as

\begin{equation}
\label{eq:heat_disc}
\theta_k = \theta_{k-1} + \tau c^{-1} f_{k-1},
\end{equation}

where $\theta_k = \theta(k \tau)$. Now, if we are given $\theta_0 = \theta(0)$, we can compute $\theta_k \; \forall k$. But what about the measurement equation given by \ref{eq:measurement_heat}? We can discretize that trivially as

\begin{equation}
\label{eq:measurement_heat_disc}
y_k = \theta_k, 
\end{equation}

where $y_k = y(k \tau)$. 

Now, consider the following problem. 

**Problem statement.** Suppose we are given the sequence $\\{y_k \\}_{k}$. Can we estimagte $\\{\theta_k \\}_k$ using the system \ref{eq:heat_disc}-\ref{eq:measurement_heat_disc}?

In a perfect world, yes of course, without any error, since then $y_k = \theta_k$. But in a world of noisy sensors and less perfect measuring equipments, we can safely assume that some error is introduced in the measurement $y_k$, i.e., instead of \ref{eq:measurement_heat_disc}, we have an equation similar to 

\begin{equation}
\label{eq:measurement_disc1}
y_k = \theta_k + v_k,
\end{equation}

where $v_k$ represents some noise. Similarly, adding energy to the object would not change its temperature "perfectly" according to \ref{eq:heat_disc}, but rather in noisy terms 

\begin{equation}
\label{eq:heat_disc1}
\theta_k = \theta_{k-1} + \tau c^{-1} f_{k-1} + w_k,
\end{equation}

where $w_k$ also represents some noise. The *problem statement* mentioned above now becomes nontrivial: given $\\{y_k\\}_k$, how can we accurately estimate $\\{\theta_k \\}_k$ using \ref{eq:measurement_disc1}-\ref{eq:heat_disc1}? 

**Example.** As an example, we consider the heating of the object over the time period $(0, 10)$ [hr]. We consider $c = 10^6$[J/m$^3$ $^\circ$ C] and $f(t) = (10 - t) 10^5$ [J / m$^3$]. We wish to estimate the temperature $\theta$ at equally spaced intervals of $1$ [hr]. 





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

We assume $v_k$ and $w_k$ to be uncorrelated Gaussian random processes with zero mean and given covariance matrices[^2]. 

**Problem statement.** We seek a solution to the following: given $\\{u_k \\}_k$ and $\\{y_k \\}_k$, how can we accurately estimate $\\{x_k \\}_k$? 

## Kalman filter

The Kalman filter is a recursive algorithm that first predicts $\hat{x_k}^{-}$ using the state equation and then uses the measurement $y_k$ to further correct $\hat{x_k}^{-}$ and produce a more accurate $\hat{x_k}^{+} = \hat{x_k} \approx x_k$. The algorithm is defined as follows. 

**1. Initialization.** Let $\hat{x_0}^{+} = x_0$ and $\Sigma_0^+ = 0$.

**2. Prediction.** Given $\hat{x_{k-1}}^{+}$, we compute $\hat{x_k}^{-}$ and $\Sigma_k^{-}$ using 

\begin{equation}
\hat{x_k}^{-} A \hat{x_{k-1}}^{+} + B u_{k-1}
\end{equation}

\begin{equation}
\Sigma_{k}^{-} = A^2 \Sigma^{+}_{k-1} + \Sigma_w
\end{equation}

**3. Correction.** In this step, we first compute the Kalman gain, $L$, as follows

\begin{equation}
L_k = \Sigma_k^- C \left(C^2 \Sigma_{k}^- + \Sigma_v \right)^{-1}.
\end{equation}

Then, we compute $\hat{x_k}^{+}$ and $\Sigma_k^{+}$ using 

\begin{equation}
\hat{x_k}^+ = \hat{x_k}^- + L_k \left(y_k - C \hat{x_k}^{-} - D u_k \right),
\end{equation}
\begin{equation}
\Sigma_k^+ = \left(I - L_k C\right) \Sigma_k^-.
\end{equation}

## References
[^1]: H.S. Carslaw, J. C. Jaeger, *Condution of Heat in Solids*, 1959, Oxford University Press.
[^2]: Chi-Tsong Chen, *Linear System Theory and Design*, 1999, Oxford University Press, 3rd Edition.
[^3]: Gregory L. Plett, *Extended Kalman filtering for battery management systems of LiPB-based HEV battery packs, Part 1, Background*, 2004, Journal of Power Sources.
