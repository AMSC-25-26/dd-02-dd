#include "schwarz_solver.hpp"


// ======================================================
// ======= TRIDIAGONAL SOLVER - Thomas algorithm ========
// ======================================================

void TridiagonalSolver::solve(const std::vector<double> &a,
                              const std::vector<double> &b,
                              const std::vector<double> &c,
                              const std::vector<double> &r,
                              std::vector<double> &x)
{
    
}


// ======================================================
// ==================== LOCAL PROBLEM ===================
// ======================================================

LocalProblem::LocalProblem(int Nnodes_global,
                           int core_start, int core_end,
                           int overlap_l,
                           double mu_, double c_,
                           double ua_, double ub_)
    : Nglob(Nnodes_global),
      core_s(core_start), core_e(core_end),
      l(overlap_l), mu(mu_), c(c_),
      ua(ua_), ub(ub_)
{
    
}


// ======================================================
// ================= SCHWARZ ITERATIONS =================
// ======================================================

SchwarzSolver::SchwarzSolver(int Nnodes_global,
                             int mpi_rank, int mpi_size,
                             int overlap_l, double mu, double c,
                             double ua, double ub,
                             int max_iter, double tol)
    : Nglob(Nnodes_global), rank(mpi_rank), size(mpi_size),
      l(overlap_l), mu(mu), c(c), ua(ua), ub(ub),
      max_iter(max_iter), tol(tol)
{

}

// Destructor: free local problem
SchwarzSolver::~SchwarzSolver() { delete local; }

void SchwarzSolver::run() {
    
}

// GATHER GLOBAL SOLUTION AND SAVE TO FILE
void SchwarzSolver::gather_and_save() {
 
}
