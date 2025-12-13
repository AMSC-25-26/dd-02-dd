# Hands-on on Domain Decomposition Technique #

This repository contains C++ implementations for solving the 1D elliptic boundary value problem

$$
\begin{aligned}
-\mu u''(x) + cu(x) = f(x), \qquad x \in (a,b),\\
\quad u(a)=u_a, \quad u(b)=u_b,
\end{aligned}
$$

using **finite difference discretization** and **domain decomposition methods**, following the course notes (`dd.pdf`).

**Advanced Methods For Scientific Computing,**

Prof. Luca Formmagia, Dr. Paolo Joseph Baioni

Polimi, Nov 2025

---

# Team 2-DD Members:

- Alessandro Russi 
- Francesca Marina Pozzi 
- Martina Rusconi 
- Micaela Perlini 
- Seyed Vahid Ghayoomie

---

# What do you have to do

For the hands-on, you should:
- implement a class to solve the problem (1). Note that you have a tridiagonal system so you can use the Thomas algorithm given in the Example. But you may choose instead to rely on the Eigen library.
- Implement the Schwarz iterator as a solver, with a class that is composed.
- implement the parallel version in MPI or OpenMP.
- have a way to shod the results (choose the graphic library or tool you prefer).

---

# Solution

## Implemented Solvers

The same discrete problem is solved using four complementary approaches:

1. **Thomas Algorithm (Sequential)**  
   - Direct solver for tridiagonal systems  
   - Serves as the baseline reference solution  

2. **Additive Schwarz Domain Decomposition**  
   - Overlapping subdomains  
   - Parallel implementation using OpenMP and MPI  
   - Iterative solver validated against the direct solution  

3. **Eigen Sparse Direct Solver**  
   - Finite difference discretization assembled into a sparse matrix  
   - Solved using Eigen’s sparse Cholesky factorization  
   - Provides an independent algebraic reference (non-FEM)  

4. **Exact Analytical Solution (Verification Benchmark)**  
   - Closed-form solutions are derived analytically for selected right-hand sides  
   - These solutions are used to validate numerical accuracy and convergence  

   **Case 1: $f(x) = \sin(\pi x)$**

   The exact solution is:
   $$
   u(x) = \frac{1}{\mu \pi^2 + c}\,\sin(\pi x)
   $$

   **Case 2: $f(x) = 1$**

   For the reaction–diffusion case ($c > 0$), the exact solution is:
   $$
   u(x) =
   \frac{1}{c}
   \left[
   1 -
   \frac{\cosh\!\left(\sqrt{\frac{c}{\mu}}\left(x-\frac{L}{2}\right)\right)}
        {\cosh\!\left(\sqrt{\frac{c}{\mu}}\frac{L}{2}\right)}
   \right]
   $$

   In the limiting Poisson case ($c = 0$), the solution reduces to:
   $$
   u(x) = \frac{x(L-x)}{2\mu}
   $$

   These analytical solutions provide an exact reference for assessing discretization errors and validating the convergence of the numerical solvers.

---

# Run

You can take a look into `2-DD.ipynb` for the solution and run it on [Google Colab 2-DD.ipynb](https://colab.research.google.com/github/AMSC-25-26/dd-02-dd/blob/main/2-DD.ipynb).
