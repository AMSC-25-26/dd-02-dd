#include <mpi.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

// ======================================================
// ======= TRIDIAGONAL SOLVER - Thomas algorithm ========
// ======================================================

class TridiagonalSolver {
public:
    static void solve(const std::vector<double> &a,
                      const std::vector<double> &b,
                      const std::vector<double> &c,
                      const std::vector<double> &r,
                      std::vector<double> &x);
};
