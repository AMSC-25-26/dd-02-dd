# Hands-on on Domain Decomposition Technique #

This repository contains C++ implementations for solving the 1D elliptic boundary value problem

$$
-\mu u''(x) + c\,u(x) = f(x), \quad x \in (a,b),\\
u(a)=u_a,\; u(b)=u_b,
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

The same discrete problem is solved using three independent approaches:

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

---

# Run

You can take a look into `2-DD.ipynb` for the solution and run it on [Google Colab 2-DD.ipynb](https://colab.research.google.com/github/AMSC-25-26/dd-02-dd/blob/main/2-DD.ipynb).
