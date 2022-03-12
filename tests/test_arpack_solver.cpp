#include "edlib/Solver/ArpackSolver.hpp"

#include <Eigen/Sparse>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include <catch2/catch.hpp>

#include <iostream>
#include <random>

constexpr uint32_t max_iter = 1000;
constexpr double tol = 1e-10;
constexpr uint32_t n_evals = 3;

template<typename MatProd>
auto eigenvectorSpectra(MatProd& prod) -> std::pair<Eigen::VectorXd, Eigen::MatrixXd>
{
    using Spectra::SortRule;
    Spectra::SymEigsSolver<MatProd> eigs(prod, n_evals, 2 * n_evals + 1);
    eigs.init();
    eigs.compute(SortRule::SmallestAlge, max_iter, tol, SortRule::SmallestAlge);
    return std::pair{eigs.eigenvalues(), eigs.eigenvectors()};
}

template<typename MatProd>
auto eigenvectorArpack(MatProd& prod) -> std::pair<Eigen::VectorXd, Eigen::MatrixXd>
{
    using namespace edlib;
    ArpackSolver solver(prod);
    solver.solve(n_evals, max_iter, tol);
    return std::pair{solver.eigenvalues(), solver.eigenvectors()};
}

TEST_CASE("Test arpack solver", "[arpack_solver]")
{
    using Catch::Matchers::WithinAbs;
    using std::abs;

    constexpr uint32_t dim = 1024;
    constexpr uint32_t n_elts = 128;

    std::mt19937 re{1337};
    std::uniform_int_distribution<uint32_t> index_dist(0, dim - 1);
    std::normal_distribution<double> elt_dist;

    Eigen::SparseMatrix<double> m(dim, dim);

    for(size_t idx = 0; idx < n_elts; ++idx)
    {
        const auto row = index_dist(re);
        const auto col = index_dist(re);
        const auto val = elt_dist(re);
        m.coeffRef(row, col) = val; // NOLINT(readability-suspicious-call-argument)
        m.coeffRef(col, row) = val; // NOLINT(readability-suspicious-call-argument)
    }
    m.makeCompressed();

    Spectra::SparseSymMatProd<double> prod(m);

    const auto [evals_spectra, evecs_spectra] = eigenvectorSpectra(prod);
    const auto [evals_arpack, evecs_arpack] = eigenvectorArpack(prod);

    REQUIRE((evals_spectra - evals_arpack).lpNorm<1>() < 1e-8);

    for(uint32_t i = 0; i < n_evals; i++)
    {
        double inner_prod = evecs_spectra.col(i).transpose() * evecs_arpack.col(i);
        REQUIRE_THAT(abs(inner_prod), WithinAbs(1.0, 1e-8));
    }
}
