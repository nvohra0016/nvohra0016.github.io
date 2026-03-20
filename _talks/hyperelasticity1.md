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
  T(u) = F(u) \left(\lamba \text{tr}(E(u)) I + 2\mu E(u) \right),
$$

where $\lambda \in \mathbb{R}$ and $\mu \in \mathbb{R}$ are Lamé parameters related to the Youngs modulus $E_Y \in \mathbb{R}$ and Poisson ratio $\nu \in \mathbb{R}$ by

$$
\label{eq:linear_fpks}
  \lambda = \frac{E_Y \nu}{(1 + \nu)(1 - 2\nu)}, \; \mu = \frac{E_Y}{2(1 + \nu)}.
$$

**Note on the choice of the St-Venant Kirchhoff material** *It can be shown that the St-Venant Kirchhoff material is not the most sound way of approximating physically observed material behaviour. Indeed, it can be shown that such a material, i.e., where the the stress relation is given by \ref{eq:st_venant_fpks}, can undergo extreme deformation to arbitrarily small volumes in a finite by expending finite energy. This, however, is physically inconsistent with natural materials. Better choices include Neo-Hookean materials, Ogden... For more information, see... Here we make use of \ref{eq:st_venant_fkps} due to its simple form.*

### 2.2.2. Linear Elasticity

In the linear elastic framework, the first Piola-Kirchhoff stress tensor is given by

$$
  T(u) = \left(\lambda \text{tr}(\epsilon(u))I + 2 \mu \epsilon(u) \right).
$$

## 2.3. Weak Formulation and Existence of Solution

We now present the weak formulation for \ref{eq:stress_eq}-\ref{eq:Dirichlet_boundary}. In this post, we consider the problem in 1D, and seek the solution $u \in H_0^1(\Omega)$. The continuous weak formulation is as follows: find $u \in H_0^1(\Omega)$ such that [^1], [Oden' 1979, Eq. 7.1] [^4],

$$
\label{eq:continuous_weak_form}
  \int_{\Omega} T(u) : \nabla \phi = \int_{\Omega} f \phi, \; \forall \phi \in H_0^1(\Omega).
$$

The existence of a solution to \ref{eq:continuos_weak_form} is established using the implicit function theorem in [Ciarlet' 1988, Theorem 6.4-1] for the St-Venant Kirchhoff material \ref{eq:st_venant_fpks}, and using Gårding operators in [Oden' 1979, Theorem 7.1] under more regularity assumptions.

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

**Existence of solution** Before we set up a numerical experiment, we first prove the existence of a solution to \ref{eq:discrete_weak_form} for appropriate forces $f$. We can rewrite \ref{eq:discrete_weak_form} as

$$
\label{eq:nonlinear_map_eq}
  \mathcal{T}(U) = L_f,
$$

where $U \in \mathbb{R}^{M-1}$ collects the entries of $u_h = \sum_{j} U_j \phi_j$, i.e., $U = [U_1 \; U_2 \; \dots \; U_{M-1}]^T$, and $\mathcal{T} : \mathbb{R}^{M-1} \rightarrow \mathbb{R}^{M-1}$ is a nonlinear map with entries

$$
  [\mathcal{T}(U)]_i = \mathcal{T}_i = \int_{\Omega} T(u_h)  \frac{d \phi_i}{dX},
$$

and the vector $L_f \in \mathbb{R}^{M-1}$ collects the entries

$$
  [L_f]_i = {L_f}_i = \int_\Omega f \phi_i.
$$

We now present a little existence result.

**Theorem 2.3.1.** Let $f \in C^0(\Omega)$. Then, for $\|f\|_\infty$ small enough, there exists a solution to \ref{eq:nonlinear_map_eq}.

*Proof* We make use of the inverse function theorem. First note that, for $U = 0$ (here we mean $0 \in \mathbb{R}^{M-1}$), we have $T(U) = 0$. Now, consider the Jacobian $J_T$ of $T$ 

$$
  [J_T(U)]_{i, j} = \frac{\partial \mathcal{T}_i}{\partial U_j}.
$$

By the chain rule, we have

$$
  [J_T(U)]_{i, j} = \int_\Omega \frac{1}{2}\left(F^2 - 1 \right) \frac{d \phi_i}{dX} \frac{d\phi_j}{dX}.
$$

## References
[^1]: Philippe G. Ciarlet, *Mathematical Elasticity: Volume 1: Three-dimensional Elasticity*, 1988, Elsevier Science Publishers.
[^2]: I. H. Shames and F. A. Cozzarelli, *Elastic and Inelastic Stress Analysis*, 1997, Taylor and Francis Group.
[^3]: J. Bonet and R. D. Wood, *Nonlinear Continuum Mechanics for Finite Element Analysis*, 1997, Cambridge University Press.
[^4]: J. Tinsley Oden, *Existence Theorems for a Class of Problems in Nonlinear Elasticity*, 1979, Journal of Mathematical Analysis and Applications.