#include "sequential_solver.hpp"
#include "schwarz_solver.hpp"
#include "thomas.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

// ======================================================
// ============== SEQUENTIAL SOLVER =====================
// ======================================================

SequentialSolver::SequentialSolver(int Nnodes, double mu_, double c_, double a_,
                                   double b_, double ua_, double ub_)
    : N(Nnodes), mu(mu_), c(c_), a(a_), b(b_), ua(ua_), ub(ub_)
{
    h = (b - a) / (N - 1);
    u.assign(N, 0.0);
    A.assign(N, 0.0);
    B.assign(N, 0.0);
    C.assign(N, 0.0);
    R.assign(N, 0.0);
}

double SequentialSolver::forcing(int i) const {
    double x = a + i * h;  // Position along the domain [a,b]
    (void)x;               // Suppress unused variable warning
    return 1.0;            // Example: constant forcing term
}

void SequentialSolver::assemble() {
    double diag_off = -mu / (h * h);
    double diag_mid = 2.0 * mu / (h * h) + c;
    
    for (int i = 0; i < N; ++i) {
        A[i] = diag_off;
        B[i] = diag_mid;
        C[i] = diag_off;
        R[i] = forcing(i);
    }
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
