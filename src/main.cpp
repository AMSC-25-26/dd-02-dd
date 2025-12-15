#include "schwarz_solver.hpp"
#include "sequential_solver.hpp"
#include <cmath>

// ======================================================================
// FORCING FUNCTIONS - Define different source terms
// ======================================================================

// 1. Constant forcing: f(x) = 1
double forcing_constant(double x) {
    (void)x;  // Suppress unused warning
    return 1.0;
}

// 2. Sinusoidal forcing: f(x) = 1 + sin(2πx)
double forcing_sin(double x) {
    return 1.0 + std::sin(2.0 * M_PI * x);
}

// 3. Parabolic forcing: f(x) = x(1-x)
double forcing_parabolic(double x) {
    return x * (1.0 - x);
}

// 4. Gaussian forcing: f(x) = exp(-50(x-0.5)²)
double forcing_gaussian(double x) {
    return std::exp(-50.0 * (x - 0.5) * (x - 0.5));
}

// 5. Piecewise forcing
double forcing_piecewise(double x) {
    if (x < 0.3) return 0.5;
    if (x < 0.7) return 2.0;
    return 0.5;
}

// 6. Polynomial forcing: f(x) = x² - x³
double forcing_polynomial(double x) {
    return x*x - x*x*x;
}

// 7. Exponential forcing: f(x) = exp(x)
double forcing_exponential(double x) {
    return std::exp(x);
}


// ======================================================================
// HELPER FUNCTION: Get forcing function name
// ======================================================================

std::string get_forcing_name(int type) {
    switch(type) {
        case 1:  return "Constant (f=1)";
        case 2:  return "Sinusoidal (f=1+sin(2πx))";
        case 3:  return "Parabolic (f=x(1-x))";
        case 4:  return "Gaussian (f=exp(-50(x-0.5)²))";
        case 5:  return "Piecewise";
        case 6:  return "Polynomial (f=x²-x³)";
        case 7:  return "Exponential (f=exp(x))";
        default: return "Unknown";
    }
}


// ======================================================================
// RELATIVE L2 ERROR COMPUTATION between two solution vectors
// ======================================================================

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


// ======================================================================
// MAIN PROGRAM
// ======================================================================

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
    int forcing_type = 1;        // forcing type selection
    
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
    if (argc >= 13) forcing_type = std::stoi(argv[12]);

    // Select forcing function based on type
    std::function<double(double)> forcing_func;

    switch (forcing_type) {
        case 1:
            forcing_func = forcing_constant;
            break;
        case 2:
            forcing_func = forcing_sin;
            break;
        case 3:
            forcing_func = forcing_parabolic;
            break;
        case 4:
            forcing_func = forcing_gaussian;
            break;
        case 5:
            forcing_func = forcing_piecewise;
            break;
        case 6:
            forcing_func = forcing_polynomial;
            break;
        case 7:
            forcing_func = forcing_exponential;
            break;
        default:
            std::cerr << "Invalid forcing type " << forcing_type 
                      << ", using constant forcing" << std::endl;
            forcing_func = forcing_constant;
            forcing_type = 1;
            break;
    }

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // ===== SEQUENTAL VERSION (only on Rank 0) =====
    std::vector<double> u_sequential;
    
    if (run_sequential && rank == 0) {
        std::cout << "\nForcing function: " << get_forcing_name(forcing_type) << std::endl;
        std::cout << "\n###############################################" << std::endl;
        std::cout << "  RUNNING SEQUENTIAL SOLVER  " << std::endl;
        std::cout << "###############################################" << std::endl;
        
        SequentialSolver seq_solver(Nnodes, mu_in, c_in, a, b, ua, ub, forcing_func);
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
                         a, b, ua, ub, max_it, tol, forcing_func);
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
