#include <mpi.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

// ======================================================
// ======= TRIDIAGONAL SOLVER - Thomas algorithm ========
// ======================================================
// 
// Input:
//   a   =   sub-diagonal, a[0] not used
//   b   =   main diagonal
//   c   =   super-diagonal, c[n-1] not used
//   r   =   right-hand side vector
//
// Output:
//   x   =   solution vector
// 
// Note: all vectors are of size n

class TridiagonalSolver {
public:
    static void solve(const std::vector<double> &a,
                      const std::vector<double> &b,
                      const std::vector<double> &c,
                      const std::vector<double> &r,
                      std::vector<double> &x);
};

// ======================================================
// ==================== LOCAL PROBLEM ===================
// ======================================================
// 
//   Nnodes_global         =  total number of nodes in global problem
//   core_start, core_end  =  global indices of local core (non-overlapping region)
//   overlap_l             =  number of overlapping nodes on each side
//   mu                    =  diffusion coefficient
//   c                     =  reaction coefficient, 
//   ua                    =  Dirichlet BC at global left boundary
//   ub                    =  Dirichlet BC at global right boundary
//   h                     =  mesh size

class LocalProblem {
public:
    LocalProblem(int Nnodes_global,
                 int core_start, int core_end,
                 int overlap_l,
                 double mu_, double c_,
                 double ua_, double ub_);

    // Forcing term f(x), right-hand side of the differential equation
    double forcing(int, double) const;

    // Assemble local tridiagonal system
    void assemble(double h);

    void apply_dirichlet(double bc_left, double bc_right);

    // Solve local tridiagonal system
    void solve_local();

    // Get internal value u[i] at global index gidx
    double value_at_global(int gidx) const;

    // Get old internal value u_old[i] at global index gidx, for communication with other processes
    double old_value_at_global(int gidx) const;

    // Get core values (non-overlapping region)
    std::vector<double> get_core_values() const;

    // Compute local squared error between current and old solution
    double local_error_sqr() const;

    // Update old solution
    void save_old();

    int ext_start() const;
    int ext_end()   const;
    int ext_length() const;
    int core_start() const;
    int core_end()   const;
    int core_length() const;

private:
    int Nglob;
    int core_s, core_e;
    int l;
    int ext_s, ext_e;
    int ext_size;
    int core_size;
    double mu, c;
    double ua, ub;

    std::vector<double> u, u_old;
    std::vector<double> A, B, C, R;
};


// ======================================================
// ================= SCHWARZ ITERATIONS =================
// ======================================================
// 
// Implements a method of overlapping domain decomposition 
// to solve a 1D problem discretized on Nglob nodes.
// Each MPI process manages a localProblem that has a "core" 
// region (its own part, non-overlapping) and an an "overlap" 
// region that extends into neighboring subdomains.
// The method iteratively solves local problems on each 
// extended subdomain, exchanging boundary values between
// neighboring processes until global convergence.
//
//   h   =   mesh size
//   bc_from_left   =   Dirichlet BC to apply at extended left boundary, got from neighbors
//   bc_from_right  =   Dirichlet BC to apply at extended right boundary, got from neighbors
//   
// Iterate until convergence or max_iter.

class SchwarzSolver {
public:
    SchwarzSolver(int Nnodes_global,
                  int mpi_rank, int mpi_size,
                  int overlap_l, double mu, double c,
                  double ua, double ub,
                  int max_iter, double tol);

    // Destructor: free local problem
    ~SchwarzSolver();

    void run();

private:
    int Nglob;
    int rank, size;
    int l;
    double mu, c;
    double ua, ub;
    int max_iter;
    double tol;

    int core_start, core_end, core_size;
    int left, right;

    LocalProblem *local;

    // Collects the local solutions and constructs 
    // u_global by averaging in the overlap points
    void gather_and_save();
};

