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
    const int n = b.size();
    x.assign(n, 0.0); 

    // Trivial cases
    if (n == 0) return;
    if (n == 1) { x[0] = r[0] / b[0]; return; }

    std::vector<double> cp(n, 0.0), dp(n, 0.0);

    cp[0] = c[0] / b[0];
    dp[0] = r[0] / b[0];

    // FORWARD ELIMINATION: eliminate sub-diagonal
    for (int i = 1; i < n; ++i) {
        // Compute modified diagonal element
        double denom = b[i] - a[i] * cp[i-1];

        // Avoid division by zero, substitute small value
        if (std::abs(denom) < 1e-18)  
            denom = (denom >= 0 ? 1e-18 : -1e-18); 
            
        // Update modified super-diagonal (except last element)
        cp[i] = (i < n-1) ? c[i] / denom : 0.0;

        // Update modified RHS
        dp[i] = (r[i] - a[i] * dp[i-1]) / denom;
    }

    // BACKWARD SUBSTITUTION: solve for x
    x[n-1] = dp[n-1];
    for (int i = n-2; i >= 0; --i)
        x[i] = dp[i] - cp[i] * x[i+1];
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
  // Extend core: every local subproblem needs l extra neighbor nodes
  // on each side (if available from neighboring subdomains).
  
  // [ext_s, ext_e] is the actual working region for this process.
  
  // The overlap allows information to propagate between subdomains:
  //   - larger l: faster convergence, more memory, more communication
  //   - smaller l: slower convergence, less memory, less communication
  
  // Minimum overlap is l = 1 (one layer of nodes)
  ext_s = std::max(0, core_s - l);            // do not go below 0
  ext_e = std::min(Nglob - 1, core_e + l);    // do not exceed Nglob-1
  ext_size = ext_e - ext_s + 1;               // number of local nodes
  core_size = core_e - core_s + 1;

  // u      =  local solution (extended)
  // u_old  =  previous solution (extended)
  //
  // Initial guess: u = 0 everywhere
  u.assign(ext_size, 0.0);
  u_old.assign(ext_size, 0.0);
    
}

// Forcing term f(x); here constant
double LocalProblem::forcing(int, double) const { return 1.0; }






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
