---
title: "An Introduction to Kalman Filters"
collection: talks
excerpt: "We present an investigation into Kalman filters and their applications in BMS"
date: 09-2025
---

## Introduction

Consider a simple 1D linear state-space system of the form 

\begin{subequations}
\begin{eqnarray}
\dot{x(t)} = A x(t),
\\
y(t) = C x(t),
\end{eqnarray}
\end{subequations}

where $x \in \mathbb{R}$ is the state variable, $y \in \mathbb{R}$ is the output, $A, C \in \mathbb{R}$ are constants. 