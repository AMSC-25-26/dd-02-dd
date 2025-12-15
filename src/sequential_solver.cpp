#include "sequential_solver.hpp"
#include "schwarz_solver.hpp"
#include "thomas.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <algorithm>

// ======================================================
// ============== SEQUENTIAL SOLVER =====================
// ======================================================

SequentialSolver::SequentialSolver(int Nnodes, double mu_, double c_, double a_,
                                   double b_, double ua_, double ub_,
                                   std::function<double(double)> forcing_func)
    : N(Nnodes), mu(mu_), c(c_), a(a_), b(b_), ua(ua_), ub(ub_), forcing(forcing_func)
{
    h = (b - a) / (N - 1);
    u.assign(N, 0.0);
    A.assign(N, 0.0);
    B.assign(N, 0.0);
    C.assign(N, 0.0);
    R.assign(N, 0.0);
}

void SequentialSolver::assemble() {
    // Initialize RHS
    R.assign(N, 0.0);

    double diag_off = -mu / (h*h);
    double diag_mid = 2.0 * mu / (h*h) + c;

    // Assign coefficients
    A.assign(N, diag_off);  // coefficients of u_{i-1}
    B.assign(N, diag_mid);  // coefficients of u_i
    C.assign(N, diag_off);  // coefficients of u_{i+1}

    // RHS: evaluate forcing at node position
    // Create a vector of indices
    std::vector<int> indices(N);
    // Fill the vector with a sequence of consecutive integers
    std::iota(indices.begin(), indices.end(), 0);

    // Apply lambda function to each element of the vector
    std::for_each(
        indices.begin(),
        indices.end(),
        [this](int i) {
            this->R[i] = this->forcing(this->a + i*this->h);
        }
    );
}

void SequentialSolver::apply_boundary_conditions() {
    // Left boundary condition
    B[0] = 1.0;
    C[0] = 0.0;
    R[0] = ua;
    
    // Right boundary condition
    B[N-1] = 1.0;
    A[N-1] = 0.0;
    R[N-1] = ub;
}

void SequentialSolver::solve() {
    assemble();
    apply_boundary_conditions();
    u = apsc::thomasSolve(B, A, C, R);
    
    std::cout << "\n================================================" << std::endl;
    std::cout << "  SEQUENTIAL SOLVER COMPLETED" << std::endl;
    std::cout << "================================================" << std::endl;
}

void SequentialSolver::save_solution(const std::string& filename) {
    std::ofstream ofs(filename);
    ofs << "x,u_s\n";
    for (int i = 0; i < N; ++i) {
        ofs << a + i * h << "," << u[i] << "\n";
    }
    ofs.close();
    std::cout << "Sequential solution saved to '" << filename << "'\n" << std::endl;
}
