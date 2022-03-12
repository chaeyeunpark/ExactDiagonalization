#include "edlib/edlib.hpp"

#include <Eigen/Dense>

#include <fstream>
#include <ios>
#include <iostream>
#include <random>

int main()
{
    using namespace edlib;
    using UINT = uint32_t;

    constexpr int N = 20;

    std::cout << "#N: " << N << std::endl;

    std::vector<double> hs;
    for(int i = 0; i <= 20; i++)
    {
        hs.emplace_back(i * 0.1);
    }
    for(auto h : hs)
    {
        std::array<Eigen::VectorXd, 2> ev;
        {
            Basis1DZ2<UINT> basis(N, 0, 1, false);
            TITFIsing<UINT> ham(basis, 1.0, h);
            const int dim = basis.getDim();

            NodeMV mv(dim, 0, dim, ham);
            auto solver = ArpackSolver(mv, dim);

            if(solver.solve(2) != ErrorType::NormalExit)
            {
                return 1;
            }
            ev[0] = solver.eigenvalues();
        }
        {
            Basis1DZ2<UINT> basis(N, 0, -1, false);
            TITFIsing<UINT> ham(basis, 1.0, h);
            const int dim = basis.getDim();

            NodeMV mv(dim, 0, dim, ham);
            auto solver = ArpackSolver(mv, dim);

            if(solver.solve(2) != ErrorType::NormalExit)
            {
                return 1;
            }
            ev[1] = solver.eigenvalues();
        }
        printf("%f\t%.10f\t%.10f\t%.10f\t%.10f\n", h, ev[0](0), ev[0](1), ev[1](0), ev[1](1));
    }

    return 0;
}
