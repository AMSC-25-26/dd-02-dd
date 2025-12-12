#ifndef SEQUENTIAL_SOLVER_HPP
#define SEQUENTIAL_SOLVER_HPP

#include <vector>
#include <string>

class SequentialSolver {
public:
    SequentialSolver(int Nnodes, double mu, double c, double ua, double ub);
    
    void solve();
    void save_solution(const std::string& filename = "sequential_solution.csv");
    
    const std::vector<double>& get_solution() const { return u; }
    
private:
    int N;                           // number of nodes
    double mu;                       // diffusion coefficient
    double c;                        // reaction coefficient
    double ua, ub;                   // boundary conditions
    double h;                        // mesh size
    
    std::vector<double> u;           // solution vector
    std::vector<double> A, B, C, R;  // tridiagonal system vectors
    
    double forcing(int i) const;
    void assemble();
    void apply_boundary_conditions();
};

#endif
