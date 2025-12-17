#include "schwarz_solver.hpp"
#include "sequential_solver.hpp"
#include <cmath>
#include <sstream>

// ======================================================================
// FORCING FUNCTIONS - Define different source terms
// ======================================================================

// 1. Constant forcing: f(x) = 1
double forcing_constant(double x) {
    (void)x;  // Suppress unused warning
    return 1.0;
}

// 2. Sinusoidal forcing: f(x) = sin(2πx)
double forcing_sin(double x) {
    return std::sin(2.0 * M_PI * x);
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
// INPUT UTILITIES
// ======================================================================

// Use of template: "this function uses a generic type T which will only be
//                   chosen when the function is called".
template <typename T>
void ask_param(const std::string& msg, T& value) {
    std::cout << msg << " [" << value << "]: ";
    std::string line;
    std::getline(std::cin, line);  // to read a line of text from input stream 

    // Default
    if (line.empty() || line == "." || line == "-") return; 

    // Convert the user’s input string into a number of type T using a stringstream, 
    // and only update the variable if the conversion succeeded
    std::stringstream ss(line);
    T tmp;
    if (ss >> tmp)
        value = tmp;
}

int ask_forcing(int def) {
    std::cout << "\nChoose forcing function:\n"
              << "  1) Constant        f(x)=1\n"
              << "  2) Sinusoidal      f(x)=sin(2πx)\n"
              << "  3) Parabolic       f(x)=x(1-x)\n"
              << "  4) Gaussian        f(x)=exp(-50(x-0.5)^2)\n"
              << "  5) Piecewise\n"
              << "  6) Polynomial      f(x)=x^2-x^3\n"
              << "  7) Exponential     f(x)=exp(x)\n"
              << "Selection [" << def << "]: ";

    std::string line;
    std::getline(std::cin, line);

    if (line.empty() || line == "." || line == "-")
        return def;

    int v = std::stoi(line);
    if (v < 1 || v > 7) {
        std::cout << "Invalid choice, using default\n";
        return def;
    }
    return v;
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
    double mu = 0.01;            // diffusion coefficient
    double c = 5.0;              // reaction coefficient
    double a = 0.0;              // left boundary of the domain
    double b = 1.0;              // right boundary of the domain
    double ua = 0.0;             // left Dirichlet BC
    double ub = 0.0;             // right Dirichlet BC
    int max_it = 2000;           // maximum number of iterations
    double tol = 1e-6;           // convergence threshold
    int forcing_type = 1;        // forcing type selection
    
    bool run_sequential = true;  // flag to run sequential solver for comparison

    if (argc >= 2) Nnodes = std::stoi(argv[1]);
    if (argc >= 3) overlap_l = std::stoi(argv[2]);
    if (argc >= 4) mu = std::stod(argv[3]);
    if (argc >= 5) c = std::stod(argv[4]);
    if (argc >= 6) a = std::stod(argv[5]);
    if (argc >= 7) b = std::stod(argv[6]);
    if (argc >= 8) ua = std::stod(argv[7]);
    if (argc >= 9) ub = std::stod(argv[8]);
    if (argc >= 10) max_it = std::stoi(argv[9]);
    if (argc >= 11) tol = std::stod(argv[10]);
    if (argc >= 12) forcing_type = std::stoi(argv[11]);
    if (argc >= 13) run_sequential = (std::stoi(argv[12]) != 0);

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        // ===== INTERACTIVE PARAMETER INPUT (only on Rank 0) =====
        if (rank == 0) {
            std::cout << "\n Press ENTER, '.' or '-' to keep default values (shown in brackets)\n\n";
    
            ask_param("Number of grid nodes", Nnodes);
            ask_param("Overlap size (in number of nodes)", overlap_l);
            ask_param("Diffusion coefficient mu", mu);
            ask_param("Reaction coefficient c", c);
            ask_param("Left domain boundary a", a);
            ask_param("Right domain boundary b", b);
            ask_param("Dirichlet left BC ua", ua);
            ask_param("Dirichlet right BC ub", ub);
            ask_param("Maximum number of iterations", max_it);
            ask_param("Tolerance", tol);
    
            int run_sequential_int = run_sequential ? 1 : 0;
            ask_param("Run sequential solver? (1=yes, 0=no)", run_sequential_int);
            run_sequential = (run_sequential_int != 0);
    
            forcing_type = ask_forcing(forcing_type);
        }
    
        // BROADCAST parameters to all ranks
        // Send the user-defined or default values from rank 0 to all other MPI processes
        MPI_Bcast(&Nnodes,            1, MPI_INT,    0, MPI_COMM_WORLD);
        MPI_Bcast(&overlap_l,         1, MPI_INT,    0, MPI_COMM_WORLD);
        MPI_Bcast(&mu,                1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&c,                 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&a,                 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&b,                 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ua,                1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ub,                1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&max_it,            1, MPI_INT,    0, MPI_COMM_WORLD);
        MPI_Bcast(&tol,               1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&forcing_type,      1, MPI_INT,    0, MPI_COMM_WORLD);
        MPI_Bcast(&run_sequential,    1, MPI_CXX_BOOL,0, MPI_COMM_WORLD);
    }

    // Forcing function selection
    std::function<double(double)> forcing;

    switch (forcing_type) {
        case 2: forcing  = forcing_sin;         break;
        case 3: forcing  = forcing_parabolic;   break;
        case 4: forcing  = forcing_gaussian;    break;
        case 5: forcing  = forcing_piecewise;   break;
        case 6: forcing  = forcing_polynomial;  break;
        case 7: forcing  = forcing_exponential; break;
        default: forcing = forcing_constant;    break;
    }
    
    // ===== SEQUENTAL VERSION (only on Rank 0) =====
    std::vector<double> u_sequential;
    
    if (run_sequential && rank == 0) {
        std::cout << "\n###############################################" << std::endl;
        std::cout << "  RUNNING SEQUENTIAL SOLVER  " << std::endl;
        std::cout << "###############################################" << std::endl;
        
        SequentialSolver seq_solver(Nnodes, mu, c, a, b, ua, ub, forcing);
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
        std::cout << "Diffusion coefficient:     " << mu << std::endl;
        std::cout << "Reaction coefficient:      " << c << std::endl;
        std::cout << "Global nodes:              " << Nnodes << std::endl;
        std::cout << "Mesh size h:               " << (b - a) / (Nnodes - 1) << std::endl;
        std::cout << "MPI processes:             " << size << std::endl;
        std::cout << "Overlap size:              " << overlap_l << " nodes" << std::endl;
        std::cout << "Tolerance:                 " << tol << std::endl;
        std::cout << "Max iterations:            " << max_it << std::endl;
        std::cout << "================================================" << std::endl;
    }

    SchwarzSolver solver(Nnodes, rank, size, overlap_l, mu, c, 
                         a, b, ua, ub, max_it, tol, forcing);
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
