**Advanced Methods For Scientific Computing**,

Prof. Luca Formaggia, Dr. Paolo Joseph Baioni

Polimi, Nov 2025

# Hands-on on Domain Decomposition Technique #

<div style="text-align: center;">
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/a9/Ddm_original_logo.png/500px-Ddm_original_logo.png" alt="dd" width="250">
</div>

This repository contains C++ implementations for solving numerically the 1D elliptic boundary value problem

$$
\begin{cases}
-\mu u''(x) + cu(x) = f(x) \qquad x \in (a,b) \\
 u(a)=u_a, \qquad u(b)=u_b
\end{cases}
$$

where:
- $\mu>0$ is the diffusion coefficient
- $c\geq0$ is the reaction coefficient
- $f(x)$ is the forcing function (source term)
- the boundary conditions are of Dirichlet type

<br />

The code uses **second order finite difference discretization** and **domain decomposition methods**, following the [course notes](docs/dd.pdf).

## Goals

1. Implement a class to solve the problem [(1)](docs/dd.pdf). 
2. Implement the Schwarz iterator as a solver, with a class that is composed.
3. Implement the parallel version in MPI or OpenMP.
4. Have a way to show the results.

## Team 2-DD Members

- **Alessandro Russi**,  &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp;11145361
- **Francesca Marina Pozzi**, &nbsp;  10837189
- **Martina Rusconi**, &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp;&nbsp;10857811
- **Micaela Perlini**, &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp;10860443
- **Seyed Vahid Ghayoomie**, &nbsp;11142478

<br />

# Code structure and overview

The project has three main components:

1. `SequentialSolver` - **Direct solver**  
    This class contains the logic needed to solve the diffusion-reaction problem in a sequential manner and is used as a reference to *validate* the parallel solver.
    - <u>Forcing</u>: various source terms are supported.
    - <u>Assembling</u>: construction of a tridiagonal matrix, stored only through its (sub-/super-)diagonal vectors.
    - <u>Solution</u>: use of *Thomas algorithm* to solve the global tridiagonal system in *linear time*.

    <br />
    
    Key methods:
    - `assemble()`, which computes the coefficients based on `mu`, `c` and grid size.
    - `apply_boundary_conditions()`, which imposes Dirichlet boundary values.
    - `solve()`, which solves the system.
    - `save_solution()`, which writes on a `.csv` files the `(x, u(x))` couples.
    
<br />

2. `LocalProblem` and `SchwarzSolver` - **Parallel MPI solver**  
    - <u>Decomposition</u>: the global domain $[a,\,b]$ is divided into `p` subdomains.
    - <u>Overlap</u>: each process extends its local domain of `l` nodes each side in order to make information exchange and convergence easier.
    - <u>Schwarz iteration</u>:
        - Each process solves its local problem in its extended region.
        - The computed boundary values falling in the overlap area are exchanged with neighboring processes using MPI.
        - The exchanged values become the new Dirichlet boundary conditions to the next iteration in the local subdomain.
    - <u>Synchronization</u>: convergence is reached when the global $\text{L}^2$ error reaches a certain tolerance.

    <br />

    The class `LocalProblem` represents the portion of global problem assigned to each MPI process.

    It distinguishes between `core_size` (= nodes owned only by a certain process) and `ext_size` (= total number of nodes on which each process works, included the overlp nodes).

    Key methods:
    - `assemble()`, which constructs the local linear system.
    - `apply_dirichlet()`, which applies boundary conditions to each loacl problem.
    - `solve_local()`, which solves the local problem using Thomas algorithm.
    - `save_old()`, which saves solution to check for convergence.
    - `local_error_sqr()`, which computes local $\text{L}^2$ error between consecutive iterations.

    <br />
    
    On the other hand, the class `SchwarzProblem` coordinates the parallel execution using MPI communication.

    In the method `run()`, it implements a loop where processes exchange boundary values and verify the global tolerance.

    The method `gather_and_save()` collects data from all the processes into Rank 0, handling averages in overlap regions to obtain a *continuous* global solution.

<br />

3. `main.cpp` - **Entry point**  
    - Handles input parameters (interactive or command line).
    - Initializes MPI.
    - Executes both sequential and parallel solver.
    - Computes relative $\text{L}^2$ error between these solutions.

## Thomas algorithm
 
A direct method to solve tridiagonal systems in **linear time**.

The <u>key method</u> found in `thomas.hpp` is `thomasSolve`, which takes as inputs four vectors of same length $N$:
- `a` = diagonal of system matrix
- `b` = sub-diagonal of system matrix
- `c` = super-diagonal of system matrix
- `f` = right-hand-side

## Configuration parameters

The program accepts the following parameters:

| Name | Symbol | Description |
| :---: | :---: | :---: |
| `Nnodes` | $N$ | Total number of nodes in the 1D grid |
| `overlap_l` | $l$ | Number of shared nodes between subdomains |
| `mu` | $\mu$ | Diffusion coefficient |
| `c` | $c$ | Reaction coefficient |
| `a`, `b` | $a,\,b$ | Domain boundaries |
| `ua`, `ub` | $u_a,\,u_b$ | Dirichlet boundary conditions |
| `max_it` | - | Maximum number of iterations |
| `tol` | - | Tolerance |
| `forcing_type` | - | Forcing type selection |
| `run_sequential` | - | Flag to run sequential solver |

**Mesh size** is given by $h=\frac{b-a}{N-1}$.

## Guide to the forcing function selection

The file `main.cpp` provides 7 different choices for the forcing function $f(x)$:

| ID | Type | Description |
| :---: | :---: | :---: |
| 1 | Constant | $f(x)=1.0$ |
| 2 | Sinusoidal | $f(x)=\sin(2\pi x)$ |
| 3 | Parabolic | $f(x)=x(1-x)$ |
| 4 | Gaussian | $f(x)=\exp(-50(x-0.5)^2)$ |
| 5 | Piecewise | Discrete values based on $x$ coordinate |
| 6 | Polynomial | $f(x)=x^2-x^3$ |
| 7 | Exponential | $f(x)=\exp(x)$ |

## Details on MPI

In Schwarz method, communication occurs in `SchwarzSolver::run()`:
1. A vector of length `L` is extracted from the overlap region of local lore.
2. Data are exchanged (simultaneous send and receive to avoid deadlocks) using `MPI_Sendrecv`.

    <u>Note</u>: processes at the boundaries of global domain are automatically handled by using `MPI_PROC_NULL`.

<br />

**Key parameter**: `overlap_l`, which is the number of shared nodes.

- `overlap_l` $\nearrow$ implies faster convergence but produces more communication overhead
- `overlap_l` $\searrow$ needs less communication but convergence is slower

The minimum overlap value is $1$.

<br />

At the end of the loop, Rank 0 gathers and saves the final solution:
1. Each process sends its extended local solution.
2. In the overlap regions, Rank 0 performs an average of the obtained local values to guarantee a *fluid transition* between subdomains.

## Output

In general, the code generates two `.csv` files:
1. `sequential_solution.csv`, which holds the result of the direct solver.
2. `parallel_solution.csv`, which holds the result of MPI parallel solver.

# Tests

Take a look into [`2-DD.ipynb`](2-DD.ipynb) to find the implementation of four complementary approaches, all to solve the same discrete problem:

1. **Thomas Algorithm (Sequential)**  
   - Direct solver for tridiagonal systems  
   - Serves as the baseline reference solution  

2. **Additive Schwarz Domain Decomposition**  
   - Overlapping subdomains  
   - Parallel implementation using OpenMP and MPI  
   - Iterative solver validated against the direct solution  

3. **Eigen Sparse Direct Solver**  
   - Finite difference discretization assembled into a sparse matrix  
   - Solved using Eigen's sparse Cholesky factorization  
   - Provides an independent algebraic reference (non-FEM)  

4. **Exact Analytical Solution (Verification Benchmark)**  

      The continuous problem reads

$$
   \begin{aligned}
   -\mu u''(x) + c u(x) = f(x), \qquad x \in (0,L), \\
   u(0) = 0, \qquad u(L) = 0 .
   \end{aligned}
$$

   Closed-form solutions are available for the following right-hand sides.

   **Case 1: $f(x) = \sin(2 \pi x)$**

$$
   \begin{aligned}
   u(x) &= \frac{1}{4 \mu \pi^2 + c} \sin(2 \pi x) .
   \end{aligned}
$$

   **Case 2: $f(x) = 1$**

   For $c > 0$:

$$
   \begin{aligned}
   u(x)
   &=
   \frac{1}{c}
   \left[
   1 -
   \frac{\cosh \left(\sqrt{\frac{c}{\mu}}\left(x-\frac{L}{2}\right)\right)}
        {\cosh \left(\sqrt{\frac{c}{\mu}}\frac{L}{2}\right)}
   \right] .
   \end{aligned}
$$

   In the limiting Poisson case $c=0$:

$$
   \begin{aligned}
   u(x) &= \frac{x(L-x)}{2\mu} .
   \end{aligned}
$$

   These exact solutions are used as verification benchmarks to assess discretization errors and validate the convergence of the numerical solvers.

<br />

The test is runnable on [Google Colab 2-DD.ipynb](https://colab.research.google.com/github/AMSC-25-26/dd-02-dd/blob/main/2-DD.ipynb), or   
[![GitHub Codespaces](https://github.com/codespaces/badge.svg)](
  https://codespaces.new/AMSC-25-26/dd-02-dd
)

# Code usage

From the project **root**, run the command
```
make -C src
```
to compile the code contained into `\src` including the header files into `\include`. This will create a new folder named `\build`, which will hold the object files and the target `schwarz_mpi`.

At this point, it is possible to run the program in different modes:
- Use
    ```
    make -C src run
    ```
    if you want to obtain **both** the sequential and the parallel MPI solutions, setting input parameters **interactively** from terminal.
- Use
    ```
    make -C src run_args
    ```
    if you want to obtain **both** the sequential and the parallel MPI solutions, but either keeping default parameters or overwriting them from **command line**.
- Use
    ```
    make -C src run_seq
    ```
    if you want to obtain the **sequential solution only**.

Each command writes `.csv` files which contain the computed solution and can be found in another new folder called `\output`.

# Visualization

To visualize the obtained results, it is suggested to use [ParaView](https://www.paraview.org) by following these instructions: 
1. Open **ParaView** on your local computer.
2. Copy the path to the `.csv` file you want to visualize.
3. Click on **File > Open** and paste the path of your file.
4. On the right side of the window, select **Line Chart View**.
5. Go to the **Properties** panel:
    - Deselect &nbsp; $\square$ **Use Index For X Axis**.
    - Choose **x** as **X Array Name**.
    - Deselect &nbsp; $\square$ **x** &nbsp; from the **Series Parameters**.
    - *Optional*: double click on the coloured circle in **Series Parameters** to choose a different color for your graph.
    - *Optional*: adjust **Bottom Axis Range** and **Left Axis Range** in such a way that your graph is visualized in a proper way.

