#include <mpi.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

// ======================================================
// ======= TRIDIAGONAL SOLVER - Thomas algorithm ========
// ======================================================
// 
// Input:
//   a   =   sub-diagonal, a[0] not used
//   b   =   main diagonal
//   c   =   super-diagonal, c[n-1] not used
//   r   =   right-hand side vector
//
// Output:
//   x   =   solution vector
// 
// Note: all vectors are of size n

class TridiagonalSolver {
public:
    static void solve(const std::vector<double> &a,
                      const std::vector<double> &b,
                      const std::vector<double> &c,
                      const std::vector<double> &r,
                      std::vector<double> &x);
};
