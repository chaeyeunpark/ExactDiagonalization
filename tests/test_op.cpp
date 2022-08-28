#include "print_vector.hpp"

#include "edlib/Op/NodeMV.hpp"
#include "edlib/Op/ParallelMV.hpp"

#include <Eigen/Sparse>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <catch2/catch_all.hpp>

#include <random>

struct SparseCol
{
    Eigen::SparseMatrix<double> m_;

    explicit SparseCol(uint32_t dim) : m_(dim, dim){};

    inline double& coeffRef(size_t row, size_t col)
    {
        return m_.coeffRef(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
    }

    [[nodiscard]] std::map<size_t, double> getCol(uint32_t col) const
    {
        using Eigen::SparseMatrix;
        std::map<size_t, double> res;

        for(SparseMatrix<double>::InnerIterator it(m_, col); it; ++it)
        {
            res.emplace(it.row(), it.value());
        }
        return res;
    }
};

TEST_CASE("Test ParallelMV", "[Op]")
{
    const size_t dim = 128;
    const size_t n_elt = 100;

    std::mt19937 re{1337};
    std::uniform_int_distribution<uint32_t> index_dist(0, dim - 1);
    std::normal_distribution<double> elt_dist;

    for(uint32_t iter = 0; iter < 100; iter++)
    {
        SparseCol m{dim};

        for(uint32_t k = 0; k < n_elt; k++)
        {
            const auto row = index_dist(re);
            const auto col = index_dist(re);
            const auto val = elt_dist(re);
            m.coeffRef(row, col) = val; // NOLINT(readability-suspicious-call-argument)
        }

        auto nmv = edlib::NodeMV(dim, 0, dim, m);
        auto pmv = edlib::ParallelMV(dim, m, 4);

        std::vector<double> vec(dim, 0.0);
        std::generate(vec.begin(), vec.end(), [&]() { return elt_dist(re); });

        std::vector<double> w1(dim, 0.0);
        std::vector<double> w2(dim, 0.0);

        nmv.perform_op(vec.data(), w1.data());
        pmv.perform_op(vec.data(), w2.data());

        std::vector<double> w_diff(dim, 0.0);
        std::transform(w1.begin(), w1.end(), w2.begin(), w_diff.begin(),
                       [](double a, double b) { return a - b; });
        INFO(w_diff);

        REQUIRE_THAT(w1, Catch::Matchers::Approx(w2));
    }
}
