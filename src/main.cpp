#include "schwarz_solver.hpp"
#include "sequential_solver.hpp"

int main(int argc, char** argv) {

    // Default parameters
    const int Nnodes = 200;            // total number of nodes in the 1D grid
    const double mu_in = 0.01;         // diffusion coefficient
    const double c_in = 5.0;           // reaction coefficient

   // ..... PARAMETERS NEEDED BY PARALLEL VERSION .....
   // -> tol, max_it, overlap_l
    
    bool run_sequential = true;  // flag to run sequential solver for comparison

    // ===== SEQUENTAL VERSION (only on Rank 0) =====
    std::vector<double> u_sequential;
    
    if (run_sequential && rank == 0) {
        std::cout << "\n###############################################" << std::endl;
        std::cout << "  RUNNING SEQUENTIAL SOLVER" << std::endl;
        std::cout << "###############################################" << std::endl;
        
        SequentialSolver seq_solver(Nnodes, mu_in, c_in, 0.0, 0.0);
        seq_solver.solve();
        seq_solver.save_solution("sequential_solution.csv");
        u_sequential = seq_solver.get_solution();
    }

    // ===== PARALLEL VERSION (all ranks) =====

    // ===== COMPARISON (only on Rank 0) =====
    
    return 0;
}
