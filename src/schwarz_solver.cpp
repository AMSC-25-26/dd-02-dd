#include "schwarz_solver.hpp"
#include "thomas.hpp"
#include <numeric>
#include <algorithm>

// ======================================================
// ==================== LOCAL PROBLEM ===================
// ======================================================

LocalProblem::LocalProblem(int Nnodes_global,
                           int core_start, int core_end,
                           int overlap_l,
                           double mu_, double c_,
                           double a_, double b_,
                           double ua_, double ub_,
                           std::function<double(double)> forcing_func)
    : Nglob(Nnodes_global),
      core_s(core_start), core_e(core_end),
      l(overlap_l), mu(mu_), c(c_), a(a_), b(b_),
      ua(ua_), ub(ub_), forcing(forcing_func)
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
            this->R[i] = this->forcing(this->a + (es + i)*h);
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

// Now u contains the solution on the extended domain [ext_s, ext_e]

// Get internal value u[i] at global index gidx
double LocalProblem::value_at_global(int gidx) const {
    if (gidx < ext_s || gidx > ext_e) return 0.0;
    return u[gidx - ext_s];
}

// Get old internal value u_old[i] at global index gidx
double LocalProblem::old_value_at_global(int gidx) const {
    if (gidx < ext_s || gidx > ext_e) return 0.0;
    return u_old[gidx - ext_s];
}

// Get core values (non-overlapping region)
std::vector<double> LocalProblem::get_core_values() const {
    std::vector<double> core(core_size);
    for (int i = 0; i < core_size; ++i) {
        int g = core_s + i;
        core[i] = u[g - ext_s];
    }
    return core;
}

// ERROR COMPUTATION
// Compute local squared error between current and old solution
double LocalProblem::local_error_sqr() const {
    double s = 0.0;
    for (int i = 0; i < ext_size; ++i) {
        double d = u[i] - u_old[i];
        s += d*d;
    }
    return s;
}

// Update old solution
void LocalProblem::save_old() { u_old = u; }

int LocalProblem::ext_start() const { return ext_s; }
int LocalProblem::ext_end()   const { return ext_e; }
int LocalProblem::ext_length() const { return ext_size; }
int LocalProblem::core_start() const { return core_s; }
int LocalProblem::core_end()   const { return core_e; }
int LocalProblem::core_length() const { return core_size; }


// ======================================================
// ================= SCHWARZ ITERATIONS =================
// ======================================================

SchwarzSolver::SchwarzSolver(int Nnodes_global,
                             int mpi_rank, int mpi_size,
                             int overlap_l, double mu, double c,
                             double a, double b,
                             double ua, double ub,
                             int max_iter, double tol,
                             std::function<double(double)> forcing_func)
    : Nglob(Nnodes_global), rank(mpi_rank), size(mpi_size),
      l(overlap_l), mu(mu), c(c), a(a), b(b), ua(ua), ub(ub),
      max_iter(max_iter), tol(tol), forcing(forcing_func)
{
  
  // DOMAIN PARTITIONING
  // Load balancing: first 'rem' processes get one extra node
  int base = Nglob / size;    // base number of nodes per rank
  int rem = Nglob % size;     // remainder nodes to distribute

  core_start = rank * base + std::min(rank, rem);  // starting global index of this process' core
  core_size  = base + (rank < rem ? 1 : 0);        // number of nodes in this process' core
  core_end   = core_start + core_size - 1;         // ending global index of this process' core (inclusive)

  // Create local problem
  local = new LocalProblem(Nglob, core_start, core_end, l, mu, c, a, b, ua, ub, forcing);

  // MPI NEIGHBOR IDENTIFICATION for excanging boundary values
  // If this process does not have left/right neighbor, set it to
  // MPI_PROC_NULL (communication with MPI_PROC_NULL is ignored by MPI)
  left  = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
  right = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

}

// Destructor: free local problem
SchwarzSolver::~SchwarzSolver() { delete local; }

void SchwarzSolver::run() {
    int Nnodes = Nglob;
    double h = (b - a) / (Nnodes - 1);

    // BOUNDARY CONDITIONS FOR LOCAL PROBLEMS
    // Will be updated at each Schwarz iteration with values from neighbors
    double bc_from_left  = ua;
    double bc_from_right = ub;

    // CONVERGENCE PARAMETERS
    double global_err = 1e9;  // initial error (as large as possible)
    int iter = 0;             // iteration counter

    while (global_err > tol && iter < max_iter) {

         // Save old solution for convergence check and communication     
        local->save_old();

        // Prepare for MPI communication of overlap values
        MPI_Status status;

        if (l > 0) {

            // Exchange a vector of lenght L
            //
            // Safety: do not send more elements than the core has 
            int L = std::min(l, core_size);

            std::vector<double> send_left_vec(L);
            std::vector<double> send_right_vec(L);

            for (int i = 0; i < L; ++i) {
                send_left_vec[i]  = local->old_value_at_global(core_start + i);
                send_right_vec[i] = local->old_value_at_global(core_end - L + 1 + i);
            }

            std::vector<double> recv_left_vec(L, ua);
            std::vector<double> recv_right_vec(L, ub);

            // MPI_Sendrecv: simultaneous send and receive (avoids deadlock)
            //
            // Send send_right_vec to the right process and receive from left in recv_left_vec
            // 
            // Safety: first Sendrecv with tag 0, second Sendrecv with tag 1
            MPI_Sendrecv(send_right_vec.data(), L, MPI_DOUBLE, right, 0,
                         recv_left_vec.data(),  L, MPI_DOUBLE, left,  0,
                         MPI_COMM_WORLD, &status);

            MPI_Sendrecv(send_left_vec.data(),  L, MPI_DOUBLE, left,  1,
                         recv_right_vec.data(), L, MPI_DOUBLE, right, 1,
                         MPI_COMM_WORLD, &status);

           
            // Select which array values to use at the boundary condition for 
            // the extended domain:
            //   - last element of the array from left neighbours
            //   - first element of the array from right neighbours
            bc_from_left  = recv_left_vec[L-1];
            bc_from_right = recv_right_vec[0];
        }
        else {
            // for l == 0, send the single conditions (last and first of the core)
            double send_left  = local->old_value_at_global(core_start);
            double send_right = local->old_value_at_global(core_end);
           
            double recv_left  = ua;
            double recv_right = ub;

            MPI_Sendrecv(&send_right, 1, MPI_DOUBLE, right, 0,
                         &recv_left,  1, MPI_DOUBLE, left,  0,
                         MPI_COMM_WORLD, &status);

            MPI_Sendrecv(&send_left,  1, MPI_DOUBLE, left,  1,
                         &recv_right, 1, MPI_DOUBLE, right, 1,
                         MPI_COMM_WORLD, &status);

            // Set to received values
            bc_from_left  = recv_left;
            bc_from_right = recv_right;
        }

        // Assemble tridiagonal system for each node
        local->assemble(h);
        local->apply_dirichlet(bc_from_left, bc_from_right);
        local->solve_local();

        // ERROR COMPUTATION
        double local_err_sq = local->local_error_sqr();
         // Combine errors from all processes using MPI_Allreduce (all processes get the result)
        MPI_Allreduce(&local_err_sq, &global_err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         // Take square root to get L2 norm
        global_err = std::sqrt(global_err);


        // Print progress (only from rank 0, every 20 iterations)
        if(rank == 0 && iter % 20 == 0){
            std::cout << "Iteration " << std::setw(4) << iter
                      << " | Error = " << std::scientific 
                      << std::setprecision(6) << global_err << std::endl;
        }

        if (!std::isfinite(global_err) || global_err > 1e300) break;

        iter++;
    }


    // Print convergence result
    if(rank == 0){
        std::cout << "================================================" << std::endl;
        if(global_err <= tol){
            std::cout << "  CONVERGED in " << iter << " iterations" << std::endl;
            std::cout << "  Final error: " << std::scientific 
                      << global_err << std::endl;
        } else {
            std::cout << "  NOT CONVERGED after " << iter << " iterations" << std::endl;
            std::cout << "  Final error: " << std::scientific 
                      << global_err << std::endl;
        }
        std::cout << "================================================" << std::endl;
    }

    gather_and_save();
}




// GATHER GLOBAL SOLUTION AND SAVE TO FILE
void SchwarzSolver::gather_and_save() {

    // Take indices from local extension
    int ext_s = local->ext_start();
    // int ext_e = local->ext_end();    (not used)
    int ext_len = local->ext_length();

    // Rank 0: reception of contributions from all processes
    if (rank == 0) {
        std::vector<double> u_global(Nglob, 0.0);
        std::vector<int> counts(Nglob, 0);

        // Rank 0 adds its own extended portion to u_global and counts
        // how many contributes each node has
        for (int i = 0; i < ext_len; ++i) {
            int g = ext_s + i;  // global index
            u_global[g] += local->value_at_global(g);
            counts[g]++;
        }

        // Receive info and the buffer value with ext_len (extended value)
        for (int p = 1; p < size; ++p) {
            MPI_Status status;
            int info[2];
            MPI_Recv(info, 2, MPI_INT, p, 100, MPI_COMM_WORLD, &status);
            int rs = info[0], rl = info[1];
            std::vector<double> buf(rl);
            MPI_Recv(buf.data(), rl, MPI_DOUBLE, p, 101, MPI_COMM_WORLD, &status);

            for (int i = 0; i < rl; ++i) {
                int g = rs + i;
                u_global[g] += buf[i];
                counts[g]++;
            }
        }

        // Average of the k contributes
        for (int i = 0; i < Nglob; ++i)
            if (counts[i] > 0) u_global[i] /= counts[i];

        // Write x and u in the global solution
        std::ofstream ofs("solution.csv");
        ofs << "x,u_p\n";
        double h = (b-a)/(Nglob-1);
        for (int i = 0; i < Nglob; ++i)
            ofs << a + i*h << "," << u_global[i] << "\n";

        ofs.close();
        std::cout << "\nSolution saved to 'solution.csv'\n";

    }
    else {
       // Other processes send their core to rank 0
        int info[2] = { ext_s, ext_len };
        MPI_Send(info, 2, MPI_INT, 0, 100, MPI_COMM_WORLD);

        std::vector<double> buf(ext_len);
        for (int i = 0; i < ext_len; ++i)
            buf[i] = local->value_at_global(ext_s + i);

        MPI_Send(buf.data(), ext_len, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD);
    }
      
}
