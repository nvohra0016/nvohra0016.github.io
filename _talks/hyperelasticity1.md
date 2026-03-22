---
title: "Into Mechanics: Investigating Hyperelasticity"
collection: talks
excerpt: "We consider hyperelasticity equations and investigate well-posedness for St-Venant Kirchhoff materials and the challenges of Newton's method as a nonlinear solver."
date: 2026-3-18
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

# 1. Introduction


# 2. Governing Equations

Let the body occupied by $\Omega$ be under some external forces, and let $\Omega'$ denote its deformed configuration (see Fig. 1 for an illustration). Since we wish to solve the system of equations to obtain the displacement field, we consider the system of equations in the reference configuration, and focus on the finite element formulation. The reader is referred to well-known references [^1] [^2] [^3] for an introduction and overview of the system of equations (the list is by no means, even vaguely, exhaustive!). We now provide a brief overview of the system of equations.

## 2.1. Kinematics

We denote the position of a particle in the reference state $\Omega$ by $X \in \Omega$. In the deformed state, we denote the position by $x = \phi (X) \in \Omega'$, where $\phi : \Omega \rightarrow \mathbb{R}^3$ denotes the deformation, and $\Omega' = \phi(\Omega)$. The reference and deformed configuration are considered in the same basis denoted by $\\{e_1, e_2, e_3 \\}$. Further, let $F(X) = \nabla \phi(X) \in \mathbb{R}^{3 \times 3}$ denote the deformation gradient, and let $u(X) = \phi(X) - I$ denote the displacement. Here and below we have $\nabla = \nabla_X$ (with respect to $X$).

We denote the right Cauchy-Green strain tensor by $C = F^T F$, and the Green-St-Venant strain tensor by $E = \frac{1}{2}\left( C - I \right)$, where $I \in \mathbb{R}^{3 \times 3}$ is the identity matrix. This gives us

$$
\label{eq:green_st_venant_strain}
  E = \frac{1}{2} \left(\nabla u + \nabla u^T + \nabla u^T \nabla u \right).
$$

Finally, the linearized strain is given by $\epsilon = \frac{1}{2}\left(\nabla u + \nabla u^T \right)$, and it can be seen that $\epsilon \approx E$ for small displacements $u$.

<div align="center">
<img src='/images/hyperelasticity1/illustration_notation.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 1. Illustration of the reference and deformed configurations.
</div>

<br>

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

In this post, we consider the simplest constitutive relation for a hyperelastic material, the St-Venant Kirchhoff relation.  For the St-Venant Kirchhoff material, the first Piola-Kirchhoff stress tensor is given by [^1]

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
  \int_{\Omega} T(u) : \nabla \psi = \int_{\Omega} f \psi, \; \forall \psi \in H_0^1(\Omega),
$$

where $:$ is the contraction operator between two tensors defined as $A : B = \sum_{i, j = 1}^{M} A_{i,j} B_{i, j}$ for any two $A, B \in \mathbb{R}^{M \times M}$. The existence of a solution to \ref{eq:continuous_weak_form} is established using the implicit function theorem in [Ciarlet' 1988, Theorem 6.4-1] for the St-Venant Kirchhoff material \ref{eq:st_venant_fpks}, and using Gårding operators in [Oden' 1979, Theorem 7.1] under more regularity assumptions.

We now present the fully discrete formulation. In this post, we focus on 1D, and in $P_1$ elements. Let $\Omega = (0, 1)$, and let $\Omega$ be divided into $M$ cells $(X_{j}, X_{j+1})$, $0 \leq j \leq M-1$, of uniform width denoted by $h = \frac{1}{M}$. Let $V_h$ is the subspace of piecewise-linear functions with basis functions given by

$$
\psi_{j}(X) = \begin{cases}
(X - X_{j-1})h^{-1}; & X \in (X_{j-1}, X_j)
\\
(X_{j+1} - X)h^{-1}; & X \in (X_j, X_{j+1})
\end{cases}, 1 \leq j \leq M-1,
$$

The fully discrete form of \ref{eq:continuous_weak_form} reads as follows: find $u_h \in V_h$ such that

$$
\label{eq:discrete_weak_form}
  \int_{\Omega} T(u_h)  \frac{d \psi_i}{dX} = \int_{\Omega} f \psi_i, \; \forall 0 \leq i \leq M.
$$

**Existence of solution.** Before we set up a numerical experiment, we first prove the existence of a solution to \ref{eq:discrete_weak_form} for appropriate forces $f$. We can rewrite \ref{eq:discrete_weak_form} as

$$
\label{eq:nonlinear_map_eq}
  \mathcal{T}(U) = L_f,
$$

where $U \in \mathbb{R}^{M-1}$ collects the entries of $u_h = \sum_{j} U_j \phi_j$, i.e., $U = [U_1 \; U_2 \; \dots \; U_{M-1}]^T$, and $\mathcal{T} : \mathbb{R}^{M-1} \rightarrow \mathbb{R}^{M-1}$ is a nonlinear map with entries

$$
  \mathcal{T}_i = \int_{\Omega} T(u_h)  \frac{d \psi_i}{dX},
$$

and the vector $L_f \in \mathbb{R}^{M-1}$ collects the entries

$$
  {L_f}_i = \int_\Omega f \psi_i.
$$

We now present a little existence result.

**Theorem 2.3.1.** Let $f \in C^0(\Omega)$. Then, for $\lVert f \rVert_\infty$ small enough, there exists a solution to \ref{eq:nonlinear_map_eq}.

*Proof.* We make use of the inverse function theorem. First note that, for $U = 0$ (here we mean $0 \in \mathbb{R}^{M-1}$), we have $\mathcal{T}(U) = 0$. Now, consider the Jacobian $\mathcal{J}$ of $\mathcal{T}$ 

$$
  {\mathcal{J}(U)}_{i, j} = \frac{\partial \mathcal{T}_i}{\partial U_j}.
$$

By definition, since $E(u_h) = \frac{1}{2}(F(u_h)^2 - 1)$, we have

$$
  \mathcal{T}_{i} = \frac{\left(\lambda + 2\mu\right)}{2}\int_\Omega \left(F(u_h)^3 - 1 \right) \frac{d\psi_i}{dX}
$$

Since $u_h = \sum_j U_j \psi_j$, by the chain rule, we have

$$
\label{eq:proof_jacobian}
  {\mathcal{J}(U)}_{i, j} = \frac{\left(\lambda + 2 \mu \right)}{2} \int_\Omega \left(3F(u_h)^2 - 1 \right) \frac{d \psi_i}{dX} \frac{d\psi_j}{dX}.
$$

For $U = 0$, we have $F(u_h) = 1$, and thus we have from \ref{eq:proof_jacobian} that $J_T(0) = \frac{\left(\lambda + 2\mu \right)}{h}\text{tri}(1, 2, 1)$ is a tri-diagonal matrix such that $\mathcal{J}_T(0)$ is symmetric positive definite. Hence $\mathcal{J}(0)$ is invertible. Thus, by the inverse function theorem [^6] $\exists$ open neighborhoods $R_1, R_2 \subset \mathbb{R}^{M-1}$ $0 \in R_1$, $0 \in R_2$, $\mathcal{T}$ is one-one on $O_1$, and

$$
  \mathcal{T}(R_1) = R_2,
$$

That is, for any $f \in R_2$, i.e., if $\lVert f \rVert_\infty$ is small enough, $\exists U_f \in R_1$ such that $\mathcal{T}(U_f) = f$. This completes the proof.

<p style="text-align: right;">&#x25A1;</p>

The above result does not prove the non-existence of solutions for any arbitrary $f$. In fact, in our numerical experiments we have obtained solutions for large $\lVert f \rVert_\infty$, however, another issue lurks with the nonlinear system above. A crucial point to consider now is how we still have not mentioned *uniqueness* for our hyperelastic system. For linear elasticity, both uniqueness and existence is well-established and follows from Korn's inequality [^1] [^5], but for hyperelastic system this is not the case. As we shall explore below, uniqueness for hyperelastic systems indeed isn't guaranteed and leads to spurious oscillations.

# 3. Numerical Experiments

## 3.1. Nonlinear Solver and Implementation Details

We now describe the details of our numerical implementation. We make use of Newton's method to solve the system \ref{eq:nonlinear_map_eq}. Starting with an initial guess $U^{(0)}$, we iterate as follows [^7]:

$$
  \mathcal{J}(U^{(m-1)}) \delta U^{(m)} = -\left(\mathcal{T}(U^{(m-1)}) - L_f \right), \nonumber
  \\
  U^{(m)} = U^{(m-1)} + \delta U^{(m)}. \nonumber
$$

We perform the Newton step till we obtain convergence of the residuals $\lVert \mathcal{T}(U^{(m-1)}) - L_f \rVert_\infty < \epsilon_{rel, tol}$ or $\lVert \delta U^{(m)} \rVert_\infty < \epsilon_{tol}$, for some prescribed tolerance $\epsilon_{abs, tol} > 0$. In our simulations, we use $\epsilon_{abs, tol} = 10^{-12}$ and $\epsilon_{rel, tol} = 10^{-14}$.

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

The implementation is done using the Python library Numpy [^9].


## 3.2. Numerical Experiment: Deformation Under Dead Load

We now consider a physical scenario with a homogenous material occupying $\Omega = (0, 1)$ [m] with $E_Y = 10^7$ [Pa] and $\nu = 0.48$ [-] to mimic the properties of rubber [^7]. We consider $f \in \\{10^6, \; 4 \times 10^7 \\}$ [N/m$^3$]. We simulate using a grid size of $h = 0.05$ [m]. We use an initial guess $U^{(0)} = 0$. The results are shown in Fig. 2.

<div align="center">
<img src='/images/hyperelasticity1/f_small_M_25.png' width='380' height='380'>
<img src='/images/hyperelasticity1/f_large_M_25.png' width='380' height='380'>
</div>

<div align = "center">
 Figure 2. Results showing the displacement due to $f = 10^6$ (left) and $f = 4 \times 10^7$ (right) forces. 
</div>

<br>

It can be observed from Fig. 2. that up to forces of O($10^6$), linear elasticity provides a good approximation and is quite close to the hyperelastic solution. For large magnitude loads, however, the difference in displacements becomes more apparent as seen in the right plot of Fig. 2. This is also consistent with the result in [Ciarlet' 1988, Theorem 6.8-1] which provides an estimate of the difference of the linear elasticity and hyperelasticity displacements up to the norm of $f$.

On the solver side, Newton's method performs well and converges within $6$ iterations when $f = 4 \times 10^6$ and within $2$ iterations when $f = 10^6$. 

**Robustness of Newton's Method.** We now investigate the robustness of Newton's method to a different initial guess, which ensures the accuracy of the solver and also helps us implement a dynamic time dependent problem later on. Since the solution is close to a parabolic profile, we choose a non-zero initial guess $U^{(0)}(X) = X(1 - X)$. The results are plotted in Fig. 3. 

<div align="center">
<img src='/images/hyperelasticity1/f_large_M_20_1.png' width='380' height='380'>
<img src='/images/hyperelasticity1/f_large_M_50_1.png' width='380' height='380'>
</div>

<div align = "center">
 Figure 3. Result showing the displacement profiles due to a non-zero initial guess for two different grid sizes: $h = 0.05$ [m] (left) and $h = 0.02$ [m] (right).
</div>

<br>

Things now take an interesting turn. The Newton solver converges, but the displacement profile is very different from what we obtained in Fig. 2. Naturally, the first instinct is to refine the grid and see what happens, but to no avail. Fig. 2 also shows the solution profile for a refined grid, which does not look promising either. Here also we have convergence of the Newton's method, albeit for around $40$ iterations this time. 

The next step is to make sure that our numerical implementation is correct, and for that reason we also implement the St-Venant Kirchhoff hyperelastic system using FEniCS [^10]. 

**Results verification using FEniCS.** We first verify our solution is using the initial guess $U^{(0)} = 0$, and then using $U^{(0)} = X(1-X)$. The results are shown in Fig. 4.

<div align="center">
<img src='/images/hyperelasticity1/f_large_M_20_fenics_comp.png' width='380' height='380'>
<img src='/images/hyperelasticity1/f_large_M_20_fenics_comp_1.png' width='380' height='380'>
</div>

<div align = "center">
 Figure 4. Results showing the displacement profiles for zero (left) and non-zero (right) initial guesses using the Author's own Numpy code implementation and using FEniCS.
</div>

<br>

For both the cases, the solution profiles are almost identical and we achieve good agreement. This gives us more confidence in our implementation, and helps us to further investigate and identify the reason for this aberrant behavior of the nonlinear solver.

# 4. Issues with Well-Posedness: Non-uniqueness of Solution

The numerical results above highlight the lack of robustness of our nonlinear solver: we have convergence, but the solution profiles are not unique. Of course, uniqueness is the first thing to look at, since the reader may have noted that till now we have only spoken of existence and no uniqueness. Indeed, we now show that we can very well have situations with our St-Venant Kirchhoff system where multiple solutions exist.

We now proceed with constructing one such scenario. Consider the Dirichlet boundary conditions 

$$
  u(0) = 0, \; u(1) = -1.5.
$$

Let $f = 0$, and let $E_Y > 0$ and $\nu \in (0, 0.5)$ be given. Consider the two solution profiles given by (see Fig. 5)

$$
u_1(X) = \begin{cases}
-X; & X \in (0, 0.5)
\\
-2X + 0.5; & X \in (0.5, 1)
\end{cases}, \;
u_2(X) = \begin{cases}
-2X; & X \in (0, 0.5)
\\
-X - 0.5; & X \in (0.5, 1)
\end{cases}.
$$

<div align="center">
<img src='/images/hyperelasticity1/non_unique_profiles.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 5. Two different displacement profiles in $H^1(\Omega)$ for boundary conditions $u(0) = 0$ and $u(1) = -1.5$.
</div>

<br>

Clearly $u_1 \in H^1(\Omega)$ and $u_2 \in H^1(\Omega)$. Now by using the definition of $F = 1 + \frac{du}{dX}$ and $E = \frac{1}{2}\left(F^2 - 1 \right)$ and \ref{eq:st_venant_fpks} we have

$$
  T(u) = \frac{(\lambda + 2\mu)}{2} \left(\frac{du}{dX}\right) \left(1 + \frac{du}{dX}\right) \left(2 + \frac{du}{dX} \right).
$$

Since $u_1$ is piecewise linear, and

$$
  \frac{du_1}{dX} = -1 \text{ on } (0, 0.5), \; \frac{du_1}{dX} = -2 \text{ on } (0.5, 1),
$$

it can be easily verified that $T(u_1) = 0$. Similar case follows for $u_2$ and it holds that $T(u_2) = 0$. Thus the variational form \ref{eq:discrete_weak_form} is satisfied for at least two different solutions, thereby providing non-uniqueness of the system. 

**Note.** *The reader may construct more such piecewise-linear displacement profiles by simply finding more roots of the polynomial*

$$
  g(\alpha) = \alpha(\alpha + 1)(\alpha + 2) - f,
$$

*for any given constant $f \in \mathbb{R}$. In fact, similar profiles can be construced for the homogeneous Dirichlet boundary conditions as well by finding $\alpha_1, \alpha_2$ which satisfy $\alpha_1 \alpha_2 < 0$ for $f > 0$, $f$ small enough, for instance.*

This helps explain the issue that we faced above, i.e., the Newton's method has performed well, and it just probably converged to a different solution *that was closer to the initial guess*. In fact, if we consider a coarse grid profile of what we have in Fig. 3 as the initial guess, then the solution converges to a similar profile for finer grids as well. 

## Further Reading and Thoughts

The above exposition provides a brief (very brief!) introduction to hyperelasticity, in particular, to nonlinear mechanics. We considered the simplest form of materials, the St-Venant Kirchhoff material, owing to its simple expression, and demonstrated the challenges of Newton's method insofar as non-uniqueness of the solution is concerned. One need not go as far as considering different intial guesses, as the author of this post has encountered issues and irregular solution profiles for $U^{(0)} = 0$ but for finer grids and for different $f$. Now we also note another issue, namely that the Jacobian may not always be invertible. Indeed, considering

$$
  {\mathcal{J}(U)}_{i, j} = \frac{\left(\lambda + 2 \mu \right)}{2} \int_\Omega \left(3F(u_h)^2 - 1 \right) \frac{d \psi_i}{dX} \frac{d\psi_j}{dX}.
$$

it can be seen when $F_h \approx \sqrt{\frac{1}{3}}$ then $\mathcal{J}$ becomes singular. This highlights the importance of a good initial guess and robustness in general. The author of this post has also investigated fixed point iteration and line search methods, but to no remarkable avail: although the fixed point iteration avoids inverting the Jacobian, it is very slow owing to its linear convergence and may require posing the system as a contraction map.

We plan to continue our investigation and consider higher dimensions and different methods and materials. That however, is the topic of a future blog post.

## References
[^1]: Philippe G. Ciarlet, *Mathematical Elasticity: Volume 1: Three-dimensional Elasticity*, 1988, Elsevier Science Publishers.
[^2]: I. H. Shames and F. A. Cozzarelli, *Elastic and Inelastic Stress Analysis*, 1997, Taylor and Francis Group.
[^3]: J. Bonet and R. D. Wood, *Nonlinear Continuum Mechanics for Finite Element Analysis*, 1997, Cambridge University Press.
[^4]: J. Tinsley Oden, *Existence Theorems for a Class of Problems in Nonlinear Elasticity*, 1979, Journal of Mathematical Analysis and Applications.
[^5]: S. C. Brenner and L. R. Scott, *The Mathematical Theory of Finite Element Methods*, 2008, Springer, 3rd Edition.
[^6]: W. Rudin, *Principles of Mathematical Analysis*, 1976, McGraw-Hill, Inc., 3rd edition.
[^7]: C. T. Kelley, *Iterative Methods for Linear and Nonlinear Equations*, 1995, Society for Industrial and Applied Mathematics.
[^8]: *Young’s Modulus of Elasticity – Values for Common Materials*, *Poisson's Ratio – Definition, Values for Materials, and Applications* Engineering Toolbox, retrieved in 2026.
[^9]: Charles R. Harris et al., *Array programming with NumPy*, 2020, Springer Science and Business Media (LLC).
[^10]: I. A. Baratta et al., *DOLFINx: The next generation FEniCS problem solving environment*, 2023.