#include "schwarz_solver.hpp"
#include "thomas.hpp"

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
  //
  // [ext_s, ext_e] is the actual working region for this process.
  //
  // The overlap allows information to propagate between subdomains:
  //   - larger l: faster convergence, more memory, more communication
  //   - smaller l: slower convergence, less memory, less communication
  //
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

// ASSEMBLE LOCAL TRIDIAGONAL SYSTEM
// h = mesh size
void LocalProblem::assemble(double h) {

    // Initialize RHS
    R.assign(ext_size, 0.0);

    double diag_off = -mu / (h*h);
    double diag_mid = 2.0 * mu / (h*h) + c;

    // Assign coefficients
    A.assign(ext_size, diag_off);   // coefficients of u_{i-1}
    B.assign(ext_size, diag_mid);   // coefficients of u_i
    C.assign(ext_size, diag_off);   // coefficients of u_{i+1}

    // RHS: evaluate forcing term at node position
    // Create a vector of local indices
    std::vector<int> local_indices(ext_size);
    // Fill the vector with a sequence of consecutive integers
    std::iota(local_indices.begin(), local_indices.end(), 0);

    // Apply lambda function to each element of the vector
    std::for_each(
        local_indices.begin(),
        local_indices.end(),
        [this, h, es = ext_s](int i) {
            int gidx = es + i;
            this->R[i] = this->forcing(gidx, h);
        }
    );
}

// APPLY BOUNDARY CONDITIONS
void LocalProblem::apply_dirichlet(double bc_left, double bc_right) {
    
    // Left boundary: 
    // IF this process owns the global left boundary, then apply Dirichlet BC
    // ELSE use the value we received from left neighbor in previous iteration
    if (ext_s == 0) {  
        B[0] = 1.0;
        C[0] = 0.0;
        R[0] = ua;
    } else {
        R[0] -= A[0] * bc_left;
        A[0] = 0.0;
    }

    // Right boundary:
    // IF this process owns the global right boundary, then apply Dirichlet BC
    // ELSE use the value we received from right neighbor in previous iteration
    if (ext_e == Nglob - 1) {
        B[ext_size-1] = 1.0;
        A[ext_size-1] = 0.0;
        R[ext_size-1] = ub;
    } else {
        R[ext_size-1] -= C[ext_size-1] * bc_right;
        C[ext_size-1] = 0.0;
    }
}

// SOLVE LOCAL TRIDIAGONAL SYSTEM
void LocalProblem::solve_local() {
    u = apsc::thomasSolve(B, A, C, R);
}

double LocalProblem::value_at_global(int gidx) const {
    if (gidx < ext_s || gidx > ext_e) return 0.0;
    return u[gidx - ext_s];
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
  // DOMAIN PARTITIONING
  // Load balancing: first 'rem' processes get one extra node
   
  int base = Nglob / size;    // base number of nodes per rank
  int rem = Nglob % size;     // remainder nodes to distribute

  core_start = rank * base + std::min(rank, rem);  // starting global index of this process' core
  core_size  = base + (rank < rem ? 1 : 0);        // number of nodes in this process' core
  core_end   = core_start + core_size - 1;         // ending global index of this process' core (inclusive)

  // Create local problem
  local = new LocalProblem(Nglob, core_start, core_end, l, mu, c, ua, ub);

  // MPI NEIGHBOR IDENTIFICATION for excanging boundary values
  // If this process does not have left/right neighbor, set it to
  // MPI_PROC_NULL (communication with MPI_PROC_NULL is ignored by MPI)
   
  left  = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
  right = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

}

// Destructor: free local problem
SchwarzSolver::~SchwarzSolver() { delete local; }

void SchwarzSolver::run() {

  // Exchange value at core boundary, not extended
  double send_left  = u[core_start - ext_start];
  double send_right = u[core_end   - ext_start];

  double recv_left  = ua;
  double recv_right = ub;

 //  First exchange: send to dest = right, receive from the src = left 
  MPI_Sendrecv(
      &send_right, 1, MPI_DOUBLE, right, 0,
      &recv_left,  1, MPI_DOUBLE, left,  0,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE
  );

 // Second exchange: send to left, receive from right
  MPI_Sendrecv(
      &send_left,  1, MPI_DOUBLE, left,  1,
      &recv_right, 1, MPI_DOUBLE, right, 1,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE
  );

  bc_left  = recv_left;
  bc_right = recv_right;
}





// GATHER GLOBAL SOLUTION AND SAVE TO FILE
void SchwarzSolver::gather_and_save() {
 
}
