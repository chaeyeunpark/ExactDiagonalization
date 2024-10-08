#include "utils.hpp"

#include "edlib/Basis/Basis1DZ2.hpp"
#include "edlib/Basis/TransformBasis.hpp"
#include "edlib/EDP/ConstructSparseMat.hpp"
#include "edlib/EDP/LocalHamiltonian.hpp"
#include "edlib/Hamiltonians/TIXXZ.hpp"
#include "edlib/Op/NodeMV.hpp"
#include "edlib/Solver/ArpackSolver.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <unsupported/Eigen/KroneckerProduct>

#include <catch2/catch_all.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>

using namespace edlib;
using Catch::Approx;

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
twoQubitOp(uint32_t N, uint32_t pos1, uint32_t pos2, // NOLINT
           const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& v1,
           const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& v2)
{
    using namespace Eigen;
    // const uint32_t dim = (1U << N);

    assert(pos1 < pos2);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> res(1, 1);
    res(0, 0) = 1.0;

    for(uint32_t i = 0; i < pos1; i++)
    {
        res = Eigen::kroneckerProduct(MatrixXd::Identity(2, 2), res).eval();
    }
    res = Eigen::kroneckerProduct(v1, res).eval();
    for(uint32_t i = pos1 + 1; i < pos2; i++)
    {
        res = Eigen::kroneckerProduct(MatrixXd::Identity(2, 2), res).eval();
    }
    res = Eigen::kroneckerProduct(v2, res).eval();
    for(uint32_t i = pos2 + 1; i < N; i++)
    {
        res = Eigen::kroneckerProduct(MatrixXd::Identity(2, 2), res).eval();
    }
    return res;
}

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
singleQubitOp(int N, int pos, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& v) // NOLINT
{
    using namespace Eigen;
    // const uint32_t dim = (1u << N);

    using MatrixT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    MatrixT res(1, 1);
    res(0, 0) = 1.0;

    for(int i = 0; i < pos; i++)
    {
        res = Eigen::kroneckerProduct(MatrixT::Identity(2, 2), res).eval();
    }

    res = Eigen::kroneckerProduct(v, res).eval();
    for(int i = pos + 1; i < N; i++)
    {
        res = Eigen::kroneckerProduct(MatrixT::Identity(2, 2), res).eval();
    }
    return res;
}

template<uint32_t N> class CompareXXZ
{
private:
    double delta_;
    Basis1DZ2<uint32_t> basis_;

    Eigen::MatrixXd hamFull_;
    double gsEnergy_;
    Eigen::VectorXd gsVec_;

public:
    constexpr static int k = (N / 2) * ((N / 2) % 2);
    constexpr static int parity = 1 - 2 * static_cast<int>((N / 2) % 2);

    explicit CompareXXZ(double delta) : delta_{delta}, basis_{N, k, parity, true}
    {
        using namespace Eigen;
        static_assert(N % 2 == 0, "N must be even");

        edp::LocalHamiltonian<double> lh(N, 2);
        for(uint32_t i = 0; i < N; i++)
        {
            lh.addTwoSiteTerm({i, (i + 1) % N}, getSXXYY() + delta_ * getSZZ());
        }
        hamFull_ = MatrixXd(edp::constructSparseMat<double>(1U << N, lh));
        SelfAdjointEigenSolver<MatrixXd> es;
        es.compute(hamFull_);

        gsEnergy_ = es.eigenvalues()[0];
        gsVec_ = es.eigenvectors().col(0);
    }

    void Test()
    {
        using namespace Eigen;
        using Catch::Approx;
        constexpr size_t max_iter = 1000;
        constexpr double tol = 1e-10;

        TIXXZ<uint32_t> ham(basis_, 1.0, delta_);
        const int dim = basis_.getDim();

        NodeMV mv(dim, 0, dim, ham);

        ArpackSolver solver(mv);
        solver.solve(2, max_iter, tol);
        const double gsEnergy1 = solver.eigenvalues()[0];

        REQUIRE(gsEnergy_ == Approx(gsEnergy1).margin(1e-4));

        const VectorXd subspaceGs = solver.eigenvectors().col(0);
        const VectorXd gsVec1 = [&]() -> VectorXd {
            auto v = toOriginalBasis(
                basis_, std::span{subspaceGs.data(), static_cast<size_t>(subspaceGs.size())});
            return Map<VectorXd>(v.data(), 1U << N);
        }();

        const double gsEnergy2
            = double(gsVec1.transpose() * hamFull_ * gsVec1) / double(gsVec1.transpose() * gsVec1);

        REQUIRE(gsEnergy1 == Approx(gsEnergy2).margin(1e-6));
        REQUIRE(std::abs(gsVec_.transpose() * gsVec1) == Approx(1.0).margin(1e-6));
    }
};

TEST_CASE("Compare GS of XXZ using LocalHamiltonian and TIBasis", "[XXZGS]")
{
    SECTION("TIBasis Z2 XXZ N=8")
    {
        CompareXXZ<8> test(1.0);
        test.Test();
    }
    SECTION("TIBasis Z2 XXZ N=10")
    {
        CompareXXZ<10> test(1.0);
        test.Test();
    }
    SECTION("TIBasis Z2 XXZ N=12")
    {
        CompareXXZ<12> test(1.0);
        test.Test();
    }
}
