---
title: "Into Mechanics: Investigating Hyperelasticity"
collection: talks
excerpt: "We consider hyperelasticity equations and investigate the robustness of Newton's method for the St Venant Kirchhoff materials."
date: 2026-3-18
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

# 1. Introduction


# 2. Governing Equations

Let the body occupied by $\Omega$ be under some external forces. Let $\Omega'$ denote the current deformed configuration. Since we wish to solve the system of equations to obtain the displacement field, we consider the system of equations in the reference configuration, and focus on the finite element formulation. The reader is referred to well-known references [^1] [^2] [^3] for... theory, and the list is by no means (even vaguely) exhaustive.

## 2.1. Kinematics

We denote the position of a particle in the reference state $\Omega$ by $X \in \Omega$. In the deformed state, we denote the position by $x = \phi (X) \in \Omega'$, where $\phi : \Omega \rightarrow \mathbb{R}^3$ denotes the deformation, and $\Omega' = \phi(\Omega)$. Further, let $F = \nabla \phi \in \mathbb{R}^{3 \times 3}$ denote the deformation gradient, and let $u(X) = \phi(X) - I$ denote the displacement. 

We denote the right Cauchy-Green strain tensor by $C = F^T F$, and the Green-St-Venant strain tensor by $E = \frac{1}{2}\left( C - I \right)$, where $I \in \mathbb{R}^{3 \times 3}$ is the identity matrix. This gives us

$$
\label{eq:green_st_venant_strain}
  E = \frac{1}{2} \left(\nabla u + \nabla u^T + \nabla u^T \nabla u \right).
$$

Finally, the linearized strain is given by $\epsilon = \frac{1}{2}\left(\nabla u + \nabla u^T \right)$, and it can be seen that $\epsilon \approx E$ for small displacements $u$.

## 2.2. Equilibrium Equations

For a pure displacement problem, under an external body force $f : \Omega \rightarrow \mathbb{R}^3$, the equation of equilibrium in the reference configuration is given by [Ciarlet' 1988, Theorem 2.6-1] [^1]

$$
\label{eq:stress_eq}
  -\nabla \cdot T(u) = f \text{ in } \Omega,
$$

where $T : \Omega \rightarrow \mathbb{R}^3$ is the second order symmetric stress tensor, and in particular the first Piola-Kirchhoff stress tensor. We complete the system \ref{eq:stress_eq} with homogeneous Dirichlet boundary conditions

$$
\label{eq:Dirichlet_boundary}
  u = 0 \text{ on } \partial \Omega.
$$

We now list the constitutive relations for the first Piola-Kirchhoff stress tensor.

### 2.2.1. Hyperelasticity: St-Venant Kirchhoff Material

In this post, we consider the simplest constitutive relation for a hyperelastic material, the St-Venant Kirchhoff relation.  For the St-Venant Kirchhoff material, the first Piola-Kirchhoff stress tensor is given by 

$$
\label{eq:st_venant_fpks}
  T(u) = F(u) \left[\lambda \text{tr}(E(u)) I + 2\mu E(u) \right],
$$

where $\lambda \in \mathbb{R}$ and $\mu \in \mathbb{R}$ are Lamé parameters related to the Youngs modulus $E_Y \in \mathbb{R}$ and Poisson ratio $\nu \in \mathbb{R}$ by

$$
\label{eq:linear_fpks}
  \lambda = \frac{E_Y \nu}{(1 + \nu)(1 - 2\nu)}, \; \mu = \frac{E_Y}{2(1 + \nu)}.
$$

**Note on the choice of the St-Venant Kirchhoff material.** *It can be shown that the St-Venant Kirchhoff material is not the most sound way of approximating physically observed material behaviour. Indeed, it can be shown that such a material, i.e., where the the stress relation is given by \ref{eq:st_venant_fpks}, can undergo extreme deformation to arbitrarily small volumes in a finite by expending finite energy. This, however, is physically inconsistent with natural materials. Better choices include Neo-Hookean materials, Ogden... For more information, see... Here we make use of \ref{eq:st_venant_fpks} due to its simple form.*

### 2.2.2. Linear Elasticity

In the linear elastic framework, the first Piola-Kirchhoff stress tensor is given by

$$
  T(u) = \left(\lambda \text{tr}(\epsilon(u))I + 2 \mu \epsilon(u) \right).
$$

## 2.3. Weak Formulation and Existence of Solution

We now present the weak formulation for \ref{eq:stress_eq}-\ref{eq:Dirichlet_boundary}. In this post, we consider the problem in 1D, and seek the solution $u \in H_0^1(\Omega)$. The continuous weak formulation is as follows: find $u \in H_0^1(\Omega)$ such that [^1], [Oden' 1979, Eq. 7.1] [^4],

$$
\label{eq:continuous_weak_form}
  \int_{\Omega} T(u) : \nabla \phi = \int_{\Omega} f \phi, \; \forall \phi \in H_0^1(\Omega),
$$

where $:$ is the contraction operator between two tensors defined as $A : B = \sum_{i, j = 1}^{M} A_{i,j} B_{i, j}$ for any two $A, B \in \mathbb{R}^{M \times M}$. The existence of a solution to \ref{eq:continuous_weak_form} is established using the implicit function theorem in [Ciarlet' 1988, Theorem 6.4-1] for the St-Venant Kirchhoff material \ref{eq:st_venant_fpks}, and using Gårding operators in [Oden' 1979, Theorem 7.1] under more regularity assumptions.

We now present the fully discrete formulation. In this post, we focus on 1D, and in $P_1$ elements. Let $\Omega = (0, 1)$, and let $\Omega$ be divided into $M$ cells $(x_{j}, x_{j+1})$, $0 \leq j \leq M-1$, of uniform width denoted by $h = \frac{1}{M}$. Let $V_h$ is the subspace of piecewise-linear functions with basis functions given by

$$
\phi_{j}(x) = \begin{cases}
(x - x_{j-1})h^{-1}; & x \in (x_{j-1}, x_j)
\\
(x_{j+1} - x)h^{-1}; & x \in (x_j, x_{j+1})
\end{cases}, 1 \leq j \leq M-1,
$$

The fully discrete form of \ref{eq:continuous_weak_form} reads as follows: find $u_h \in V_h$ such that

$$
\label{eq:discrete_weak_form}
  \int_{\Omega} T(u_h)  \frac{d \phi_i}{dX} = \int_{\Omega} f \phi_i, \; \forall 0 \leq i \leq M.
$$

**Existence of solution.** Before we set up a numerical experiment, we first prove the existence of a solution to \ref{eq:discrete_weak_form} for appropriate forces $f$. We can rewrite \ref{eq:discrete_weak_form} as

$$
\label{eq:nonlinear_map_eq}
  \mathcal{T}(U) = L_f,
$$

where $U \in \mathbb{R}^{M-1}$ collects the entries of $u_h = \sum_{j} U_j \phi_j$, i.e., $U = [U_1 \; U_2 \; \dots \; U_{M-1}]^T$, and $\mathcal{T} : \mathbb{R}^{M-1} \rightarrow \mathbb{R}^{M-1}$ is a nonlinear map with entries

$$
  \mathcal{T}_i = \int_{\Omega} T(u_h)  \frac{d \phi_i}{dX},
$$

and the vector $L_f \in \mathbb{R}^{M-1}$ collects the entries

$$
  {L_f}_i = \int_\Omega f \phi_i.
$$

We now present a little existence result.

**Theorem 2.3.1.** Let $f \in C^0(\Omega)$. Then, for $\lVert f \rVert_\infty$ small enough, there exists a solution to \ref{eq:nonlinear_map_eq}.

*Proof.* We make use of the inverse function theorem. First note that, for $U = 0$ (here we mean $0 \in \mathbb{R}^{M-1}$), we have $\mathcal{T}(U) = 0$. Now, consider the Jacobian $\mathcal{J}$ of $\mathcal{T}$ 

$$
  {\mathcal{J}(U)}_{i, j} = \frac{\partial \mathcal{T}_i}{\partial U_j}.
$$

By definition, since $E(u_h) = \frac{1}{2}(F(u_h)^2 - 1)$, we have

$$
  \mathcal{T}_{i} = \frac{\left(\lambda + 2\mu\right)}{2}\int_\Omega \left(F(u_h)^3 - 1 \right) \frac{d\phi_i}{dX}
$$

Since $u_h = \sum_j U_j \phi_j$, by the chain rule, we have

$$
\label{eq:proof_jacobian}
  {\mathcal{J}(U)}_{i, j} = \frac{\left(\lambda + 2 \mu \right)}{2} \int_\Omega \left(3F(u_h)^2 - 1 \right) \frac{d \phi_i}{dX} \frac{d\phi_j}{dX}.
$$

For $U = 0$, we have $F(u_h) = 1$, and thus we have from \ref{eq:proof_jacobian} that $J_T(0) = \frac{\left(\lambda + 2\mu \right)}{h}\text{tri}(1, 2, 1)$ is a tri-diagonal matrix such that $\mathcal{J}_T(0)$ is symmetric positive definite. Hence $\mathcal{J}(0)$ is invertible. Thus, by the inverse function theorem [^5] $\exists$ open neighborhoods $R_1, R_2 \subset \mathbb{R}^{M-1}$ $0 \in R_1$, $0 \in R_2$, $\mathcal{T}$ is one-one on $O_1$, and

$$
  \mathcal{T}(R_1) = R_2,
$$

That is, for any $f \in R_2$, i.e., if $\lVert f \rVert_\infty$ is small enough, $\exists U_f \in R_1$ such that $\mathcal{T}(U_f) = f$. This completes the proof.

<p style="text-align: right;">&#x25A1;</p>

The above result does not prove the non-existence of solutions for any arbitrary $f$. In fact, in our numerical experiments we have obtained solutions for large $\lVert f \rVert_\infty$, however, another issue lurks with the nonlinear system above. A crucial point to consider now is how we still have not mentioned *uniqueness* for our hyperelastic system. For linear elasticity, both uniqueness and existence is well-established and follows from Korn's inequality [^1], but for hyperelastic system this is not the case. As we shall explore below, uniqueness for hyperelastic systems indeed isn't guaranteed and leads to spurious oscillations.

# 3. Numerical Experiments

## 3.1. Nonlinear Solver and Implementation Details

We now describe the details of our numerical implementation. We make use of Newton's method to solve the system \ref{eq:nonlinear_map_eq}. Starting with an initial guess $U^{(0)}$, we perform the step [^6]

$$
  \mathcal{J}(U^{(m-1)}) \delta U^{(m)} = -\left(\mathcal{T}(U^{(m-1)}) - L_f \right), \nonumber
  \\
  U^{(m)} = U^{(m-1)} + \delta U^{(m)}. \nonumber
$$

We perform the Newton step till we obtain convergence of the residuals $\lVert \mathcal{T}(U^{(m-1)}) - L_f \rVert_\infty < \epsilon_{tol}$ or $\lVert \delta U^{(m)} \rVert_\infty < \epsilon_{tol}$, for some prescribed tolerance $\epsilon_{tol} > 0$.

Since we are using $P_1$ elements for the displacement $u$, this means that $F$ is a piecewise-constant on each grid cell. This makes numerical integration straightforward when computing the Jacobians of $\mathcal{T}$. Indeed, from \ref{eq:proof_jacobian} we have

$$
  \mathcal{J}_{i, i} = \frac{\left(\lambda + 2\mu \right)}{2h}\left( 3F_{i}^2 + 3F_{i+1}^2 - 2  \right), \; F_i = \frac{\left(U_{i+1} - U_i \right)}{2}, \; F_{i+1} = \frac{\left(U_{i+2} - U_{i+1} \right)}{2}, \nonumber
$$

$$
  \mathcal{J}_{i, i+1} = -\frac{\left(\lambda + 2\mu \right)}{2h}\left(3F_{i+1}^2 - 1  \right), \; F_i = \frac{\left(U_{i+1} - U_i \right)}{2}, \nonumber
$$

$$
  \mathcal{J}_{i, i-1} = -\frac{\left(\lambda + 2\mu \right)}{2h}\left(3F_{i}^2 - 1  \right), \; F_i = \frac{\left(U_{i+1} - U_i \right)}{2}. \nonumber
$$

## 3.2. Numerical Experiment: Deformation Under Dead Load

We now consider a physical scenario with a homogenous material occupying $\Omega = (0, 1)$ [m] with $E_Y = 10^6$ [Pa] and $\nu = 0.45$ [-]. We consider $f = 100$. We simulate using a grid size of $h = 0.05$ [m].

## References
[^1]: Philippe G. Ciarlet, *Mathematical Elasticity: Volume 1: Three-dimensional Elasticity*, 1988, Elsevier Science Publishers.
[^2]: I. H. Shames and F. A. Cozzarelli, *Elastic and Inelastic Stress Analysis*, 1997, Taylor and Francis Group.
[^3]: J. Bonet and R. D. Wood, *Nonlinear Continuum Mechanics for Finite Element Analysis*, 1997, Cambridge University Press.
[^4]: J. Tinsley Oden, *Existence Theorems for a Class of Problems in Nonlinear Elasticity*, 1979, Journal of Mathematical Analysis and Applications.
[^5]: W. Rudin, *Principles of Mathematical Analysis*, 1976, McGraw-Hill, Inc., 3rd edition.
[^6]: C. T. Kelley, *Iterative Methods for Linear and Nonlinear Equations*, 1995, Society for Industrial and Applied Mathematics.