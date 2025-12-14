#include "schwarz_solver.hpp"
#include "sequential_solver.hpp"

// Function to compute relative L2 error between the two solutions
double compute_relative_error(const std::vector<double>& u1, 
                              const std::vector<double>& u2) {
    if (u1.size() != u2.size()) return -1.0;
    
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < u1.size(); ++i) {
        double diff = u1[i] - u2[i];
        num += diff * diff;
        den += u2[i] * u2[i];
    }
    return std::sqrt(num / den);
}

int main(int argc, char** argv) {

    // Default parameters
    int Nnodes = 200;            // total number of nodes in the 1D grid
    int overlap_l = 4;           // number of shared nodes between domains
    double mu_in = 0.01;         // diffusion coefficient
    double c_in = 5.0;           // reaction coefficient
    double a = 0.0;              // left boundary of the domain
    double b = 1.0;              // right boundary of the domain
    double ua = 0.0;             // left Dirichlet BC
    double ub = 0.0;             // right Dirichlet BC
    int max_it = 2000;           // maximum number of iterations
    double tol = 1e-6;           // convergence threshold
    
    bool run_sequential = true;  // flag to run sequential solver for comparison

    // Take input values if user passes them as arguments
    if (argc >= 2) Nnodes = std::stoi(argv[1]);
    if (argc >= 3) overlap_l = std::stoi(argv[2]);
    if (argc >= 4) mu_in = std::stod(argv[3]);
    if (argc >= 5) c_in = std::stod(argv[4]);
    if (argc >= 6) a = std::stoi(argv[5]);
    if (argc >= 7) b = std::stoi(argv[6]);
    if (argc >= 8) ua = std::stoi(argv[7]);
    if (argc >= 9) ub = std::stoi(argv[8]);
    if (argc >= 10) max_it = std::stoi(argv[9]);
    if (argc >= 11) tol = std::stod(argv[10]);
    if (argc >= 12) run_sequential = (std::stoi(argv[11]) != 0);

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // ===== SEQUENTAL VERSION (only on Rank 0) =====
    std::vector<double> u_sequential;
    
    if (run_sequential && rank == 0) {
        std::cout << "\n###############################################" << std::endl;
        std::cout << "  RUNNING SEQUENTIAL SOLVER  " << std::endl;
        std::cout << "###############################################" << std::endl;
        
        SequentialSolver seq_solver(Nnodes, mu_in, c_in, a, b, ua, ub);
        seq_solver.solve();
        seq_solver.save_solution("sequential_solution.csv");
        u_sequential = seq_solver.get_solution();
    }

    // Synchronization before parallel run
    MPI_Barrier(MPI_COMM_WORLD);

    // ===== PARALLEL VERSION (all ranks) =====
    if(rank == 0) {
        std::cout << "\n###############################################" << std::endl;
        std::cout << "  RUNNING PARALLEL SCHWARZ SOLVER" << std::endl;
        std::cout << "###############################################" << std::endl;
        std::cout << "\n================================================" << std::endl;
        std::cout << "  SCHWARZ METHOD CONFIGURATION\n";
        std::cout << "================================================" << std::endl;
        std::cout << "Diffusion coefficient:     " << mu_in << std::endl;
        std::cout << "Reaction coefficient:      " << c_in << std::endl;
        std::cout << "Global nodes:              " << Nnodes << std::endl;
        std::cout << "Mesh size h:               " << (b - a) / (Nnodes - 1) << std::endl;
        std::cout << "MPI processes:             " << size << std::endl;
        std::cout << "Overlap size:              " << overlap_l << " nodes" << std::endl;
        std::cout << "Tolerance:                 " << tol << std::endl;
        std::cout << "Max iterations:            " << max_it << std::endl;
        std::cout << "================================================" << std::endl;
    }

    SchwarzSolver solver(Nnodes, rank, size, overlap_l, mu_in, c_in, 
                         a, b, ua, ub, max_it, tol);
    solver.run();

    // ===== COMPARISON (only on Rank 0) =====
    if (rank == 0) {
        std::cout << "\n###############################################" << std::endl;
        std::cout << "  PERFORMANCE COMPARISON  " << std::endl;
        std::cout << "###############################################" << std::endl;
        
        if (run_sequential) {
            std::vector<double> u_parallel(Nnodes);
            std::ifstream parallel_file("parallel_solution.csv");
            std::string header;
            std::getline(parallel_file, header);
            
            for (int i = 0; i < Nnodes; ++i) {
                double x, u_p;
                char comma;
                parallel_file >> x >> comma >> u_p;
                u_parallel[i] = u_p;
            }
            parallel_file.close();
            
            double rel_error = compute_relative_error(u_parallel, u_sequential);
            std::cout << "Relative L2 error:         " << std::scientific 
                      << std::setprecision(6) << rel_error << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
