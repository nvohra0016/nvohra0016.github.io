---
title: "Accelerating Convergence of Newton's Method for Degenerate Functions"
collection: talks
excerpt: "We explore the application of Newton's method when solving the Butler-Volmer equation, and ways to reduce the number of iterations taken by the solver."
date: 2025-10-19
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

# 1. Introduction to Newton's method
Newton's method is a perpetual workhorse that has proven to be one of the most successful and versatile nonlinear solvers. The framework has been applied to cover a large range of nonlinearity, from s to piecewise-smooth (semismooth) functions [^1] [^2]. 

## 1.1 Algorithm 
Let $F \in C^\infty(\mathbb{R}), \; F : \mathbb{R} \rightarrow \mathbb{R}$ be a given smooth function, and we wish to solve $F(x) = 0$. The Newton's method generates a sequence $\\{ x^{(m)} \\}$ iteratively: given $x^{(m-1)}$, we obtain $x^{(m)}$ as

\begin{equation}
\label{eq:Newton_method1}
    R^{(m-1)} = F \left(x^{(m-1)} \right),
\end{equation}

\begin{equation}
\label{eq:Newton_method2}
    \delta^{(m-1)} = - {J^{(m-1)}}^{-1} R^{(m-1)},
\end{equation}

\begin{equation}
\label{eq:Newton_method3}
    x^{(m)} = x^{(m-1)} + \delta^{(m-1)},
\end{equation}

where $J^{(m-1)} = F'\left( x^{(m-1)} \right)$ is the Frechet derivative (the Jacobian of $F$), and $x^{(0)} = x_0$ is the inital guess which we are given. We iterate till a prescribed tolerance $\epsilon_{tol} > 0$ is reached by the absolute of the residual $R^{(m)}$.

**Note on the choice of Jacobian.** *In the semismooth framework, the Jacobian $J \in \partial_B F(x)$ is the Clarke's generalized Jacobian which is computed using the B-subdifferential*

\begin{equation}
\partial_B F(x) = \\{ J_F \in \mathbb{R} \; | \; \exists \{x_k\} \in D_F, \; x_k \rightarrow x, \; \left(F \right)'(x_k) \rightarrow J_F \\},
\end{equation}

 *where for each $x \in D_F \subset \mathbb{R}$, the Fre´chet derivative $F'(x)$ exists [^2]*.
 

## 1.2 Convergence of Newton's method
Under the assumptions that $F$ is Lipschitz with a bounded derivative, convergence of the algorithm \ref{eq:Newton_method1}-\ref{eq:Newton_method3} has been established for an appropriate initial guess $x^{(0)}$[^1]. Here we present a short proof adapted from[^2] for smooth functions.

**Theorem** Let $F \in C^\infty (\mathbb{R})$ and let $x_{sol}$ be the unique solution $F(x_{sol}) = 0$. Suppose $J_x = F'(x) \neq 0, \; \forall x \in \mathbb{R}$. Then...

*Proof.*  Let $\delta > 0$ and let

\begin{equation}
    \left|J_x^{-1} \right| \leq M, \; \left|F''(x) \right| \leq C, \; \forall x \in (x_{sol} - \delta, x_{sol} + \delta).
\end{equation}

Since $F$ is smooth, by Taylor's theorem $\exists \beta \in (x_{sol}-\delta, x_{sol} + \delta)$ such that[^7]

\begin{equation}
\label{eq:proof1}
    \left|F(x_{sol}) - F(x_{sol} - \delta) - J_{(x_{sol} - \delta)} \delta \right| \leq \frac{|F''(\beta) \delta^2|}{2} \leq C \left| \delta \right|^2.
\end{equation}

Let $v^{(m-1)} = x^{(m-1)} - x_{sol}$. Thus $v^{(m)} - v^{(m-1)} = \delta^{(m-1)}$. Then, from \ref{eq:Newton_method2} we have

\begin{equation}
\label{eq:proof2}
 -J^{(m-1)}\left(v^{(m)} - v^{(m-1)} \right) = F \left(x_{sol} + v^{(m-1)} \right) - F(x_{sol}),
 \end{equation}

since $F(x_{sol}) = 0$. Rewriting \ref{eq:proof2} we get

\begin{equation}
\label{eq:proof3}
    -J^{(m-1)} v^{(m)} = F\left(x_{sol} + v^{(m-1)} \right) - F(x_{sol}) - J^{(m-1)} v^{(m-1)}.
\end{equation}

Finally, using the estimate from \ref{eq:proof1}, we have

\begin{equation}
\label{eq:proof4}
    \left| v^{(m)} \right| = \left| \left( {J^{(m)}} \right)^{-1} J^{(m)} v^{(m)} \right| \leq M C \left|v^{(m-1)} \right|^2.
\end{equation}

Or, by induction, we have from \ref{eq:proof4}

\begin{equation}
    \left| v^{(m)} \right| \leq \left(MC \right)^{2^m - 1} \left|v^{(0)} \right|^{2^m}.
\end{equation}

# 2. Butler-Volmer equation

The Butler-Volmer equation is used to model the electrochemical reaction kinetics taking place in electrochemical batteries. In particular, they are used to model the lithium ion exchange between the electrode and electrolyte at the microscopic and macroscopic scale[^3]. The equation describes the relationship between the current density $j$ [A/m$^2$] and the electric potential $\phi$ [V] as[^4] 

\begin{equation}
\label{eq:BV}
    j(\phi) = j_0 \left(e^{\frac{\alpha_a z F}{RT} \phi} -  e^{-\frac{\alpha_c z F}{RT} \phi}\right),
\end{equation}

where $j_0$ [A/m] is the exchange current density, $\alpha_a$ and $\alpha_c$ [-] are the cathodic and anodic charge transfer coefficients such that $\alpha_a + \alpha_c = 1$, $F \approx 9.648 \times 10^4$ [C/mol] is Faraday's constant, $R \approx 8.314$ [J/K mol] is the gas constant, $z$ is the number of electrons involved in the electrode reaction (for ex. $z = 1$), and $T$ [K] is the temperature. Fig. 1 shows a plot of the current density as a function of the voltage for commonly chosen physical parameters and $j_0 = 1$ [A/m].

<div align="center">
<img src='/images/Newtons_method_images/BV.png' width='380' height='380'>
<img src='/images/Newtons_method_images/BV1.png' width='380' height='380'>
</div>

<div align = "center">
 Figure 1. Plot showing the current density as a function of the potential using the Butler-Volmer equation. Here $\alpha_c = \alpha_a = 0.5$. Left: for $\phi \in [-0.1, 0.1]$. Right: for $\phi \in [-0.5, 0.5]$.
</div>

<br>

## 2.1 Example: solving the Butler-Volmer equation
We now consider solving the 1D equation

\begin{equation}
\label{eq:example_BV}
    j(\phi) + \sigma \phi = c,
\end{equation}

where $\sigma$ [S/m] is the conductivity, and the current density $j$ is described by \ref{eq:BV}. For the purposes of this example, we consider $\sigma = 10^6$ to represent the physical values of most materials[^5]. The value of $c$ in the RHS of \ref{eq:example_BV} will be chosen so as to make the solution $\phi_{sol}$ of \ref{eq:example_BV} closer towards $1$ (which increases the gradient of $j$!). The tolerance is set to $\epsilon_{tol} = 10^{-6}$. 

We test with multiple values $c \in \\{10^5, 5 \times 10^5, 10^6 \\}$. We consider an initial guess of $\phi^{(0)} = 0$. The absolute value of the residuals $R = j(\phi) +\sigma \phi - c$ is plotted in Fig. 2 (left) below. The solution for $c = 10^6$ is $\phi_{sol} = 0.6549109190$ and $j_{sol} = 3.450890810218518 \times 10^5$.

<div align="center">
<img src='/images/Newtons_method_images/alpha_res_smooth.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 2. Plot showing the residuals for different values of $c$ when $\alpha_a = 0.5$.
</div>

<br>

It can be observed that increasing the value of $c$ increases the number of iterations from $3$ to $5$ to $13$. This happens because increasing $c$ increases the potential $\phi$, which further increases the gradient $J$, leading to nondegenerate behaviour.  

**Note.** *Equations like \ref{eq:example_BV} are encountered as boundary conditions after spatially discretizing the governing equations using, for example, the finite element or finite volume method. Presently, we do not discuss any numerical discretizations, but will revisit that in a future blog post!*

## 2.2 Example: non-convergence for large charge transfer coefficients

We now consider $\alpha_a = 0.85$ in \ref{eq:BV} in order to test the robustness of Newton's method; see Fig. 3 below for a plot of the current density. It can be observed that increasing $\alpha_a$ increases the gradients near $\phi = 1$.

<div align="center">
<img src='/images/Newtons_method_images/BV_alpha85.png' width='380' height='380'>
<img src='/images/Newtons_method_images/BV1_alpha85.png' width='380' height='380'>
</div>

<div align = "center">
 Figure 3. Plot showing the current density as a function of the potential using the Butler-Volmer equation. Here $\alpha_a = 0.85, \; \alpha_c = 0.15$. Left: for $\phi \in [-0.1, 0.1]$. Right: for $\phi \in [-0.5, 0.5]$.
</div>

<br>

We now return to solving \ref{eq:example_BV} using $c = 10^6$. In this case, the number of iterations taken by the solver increases to $26$. The reported solution is $\phi_{sol} = 0.4018586997$ and $j_{sol} = 5.981413003436504 \times 10^5$. The convergence can be improved by choosing a different initial guess. For example, for $\phi^{(0)} = 0.5$, the number of iterations taken drops to $9$ and the reported solution is $\phi_{sol} = 0.4018586997$ and $j_{sol} = 5.981413003430553 \times 10^5$; see Fig. 4 for a comparison of the residuals.

<div align="center">
<img src='/images/Newtons_method_images/alpha_res_degenerate.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 4. Plot showing the residuals for different values of the initial guess $\phi^{(0)}$ when $\alpha_a = 0.85$. Note that the convergence improves slightly when the initial guess is chosen such that it is close to the solution.
</div>


## 2.2 Improving convergence of Newton's method: primary variable switch

We now seek to improve the convergence of the Newtons method without having to fine tune our initial guess too much. To this end, we introduce a new formulation of \ref{eq:example_BV} with $j$ as the primary variable as follows: find a solution $j_{sol}$ to

\begin{equation}
    \label{eq:example_BV_j}
    j + \sigma \phi(j) = c,    
\end{equation}

where $\phi(j)$ represents the inverse of $j(\phi)$, i.e., we switch the primary variable from $\phi$ to $j$. Fig. 5 shows a plot of the potential $\phi$ as a function of the current density $j$ using \ref{eq:BV} for $\alpha_a = \alpha_c = 0.5$.

<div align="center">
<img src='/images/Newtons_method_images/BV_beta.png' width='380' height='380'>
<img src='/images/Newtons_method_images/BV_beta1.png' width='380' height='380'>
</div>

<div align = "center">
 Figure 5. Plot showing the potential as a function of the current density using the Butler-Volmer equation. Here $\alpha_c = \alpha_a = 0.5$. Left: for $j \in [-6, 6]$. Right: for $j \in [-15000, 15000]$.
</div>

<br>

We now discuss the motivation for changing our primary variable from $\phi$ to $j$. As observed earlier, when $\alpha_a = 0.85$, the Butler-Volmer equation leads to a high gradient as $\phi$ increases.
In particular, for the case of $\alpha_a = 0.85$, the function $\phi(j)$ exhibits concave behaviour with bounded derivatives near $j = \phi(1)$, as opposed to $j(\phi)$ near $\phi = 1$ as seen in Fig. 3. This is made clear in Fig. 6 below.

<div align="center">
<img src='/images/Newtons_method_images/BV_degenerate_beta.png' width='380' height='380'>
<img src='/images/Newtons_method_images/BV_degenerate_beta1.png' width='380' height='380'>
</div>

<div align = "center">
 Figure 6. Plot showing the potential $\phi(j)$ as a function of the current density near $j = \phi(1)$. Here $\alpha_c = \alpha_a = 0.5$. Left: for $j \in [-1, 30]$. Right: for $j \in [-1, 2 \times 10^8]$. Note that the derivative $\phi'$ remains bounded and small as $j$ becomes large (as $\phi$ approaches 1).
</div>

<br>

We now solve \ref{eq:example_BV_j} with initial guess $j^{(0)} = \phi(0)$ and $j^{(0)} = \phi(0.5)$ to check if the change in the primary variable improves the convergence of the Newton's method. The reported solution in both cases are $\phi_{sol} = 0.4018586997$,$j_{sol} = 5.981413003440497 \times 10^5$ and $\phi_{sol} = 0.4018586997$,$j_{sol} = 5.981413003440507 \times 10^5$, respectively. Fig. 7 shows a plot comparing the residuals for different initial guesses comparing the two primary variable approaches.

<div align="center">
<img src='/images/Newtons_method_images/beta_alpha_res_degenerate.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 7. Plot showing the residuals computed using primary variable (abbreviated to Pri. var. in legend) $\phi$ and $j$ when $\alpha_a = 0.85$. Note the improvement in the convergence when $j$ is used as a primary variable.
</div>

<br>

Indeed, it can be observed that for the case when $j^{(0)} = 0$, the number of iterations have been reduced to from $26$ to $8$. When $j^{(0)} = \phi(0.5)$, the number of iterations taken are also $8$. This example demonstrates the robustness of using $j$ as the primary variable when solving the Butler-Volmer equation.

**Note on computing the inverse function $\phi(j)$.** *The reader may have guessed that the inverse of \ref{eq:BV} does not have a simple analytical expression. Indeed, for our implementation, we have made use of the Ridder method as part of the SciPy library[^6] to aid in computing the inverse $\phi(j)$, i.e., given $j_0 \in \mathbb{R}$, we solve for $\phi_0$ in $j_0 - \phi(\phi_0) = 0$ where $\phi$ is given by \ref{eq:BV}. It may be noted that Newton's method may be used instead of Ridder method as well without affecting the results presented here.*


## References
[^1]: C.T. Kelley, *Iterative Methods for Linear and Nonlinear Equations*, 1995, Society for Industrial and Applied Mathematics.
[^2]: Michael Ulbrich, *Semismooth Newton Methods for Variational Inequalities and Constrained Optimization Problems in Function Spaces*, 2011, Mathematical Optimization Society and the Society for Industrial and Applied Mathematics.
[^3]: F. Brosa et al., *A continuum of physics-based lithium-ion battery models reviewed*, 2022, Progress in Energy, 4.
[^4]: Gregory L. Plett, *Battery Management Systems, Volume 1: Battery Modeling Battery Modeling*, 2015, Artech House Publishers.
[^5]: *Electrical Conductivity of some Common Materials*, Engineering Toolbox, retreived in 2025.
[^6]: P. Virtanen et al., *SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python*, 2020, Nature Methods, 17.
[^7]: W. Rudin, *Principles of Mathematical Analysis*, 1976, McGraw-Hill, Inc., 3rd edition.