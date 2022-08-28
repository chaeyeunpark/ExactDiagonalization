#include "edlib/Basis/Basis1D.hpp"
#include "edlib/Basis/Basis1DZ2.hpp"
#include "edlib/EDP/ConstructSparseMat.hpp"
#include "edlib/Hamiltonians/TIXXZ.hpp"

#include "XXZ.hpp"
#include "utils.hpp"

#include <Eigen/Dense>
#include <catch2/catch_all.hpp>

#include <iostream>
#include <unordered_set>

using namespace Eigen;
using namespace edlib;

constexpr uint32_t max_n = 12;

template<typename Basis> VectorXd translate(Basis&& basis, const VectorXd& r)
{
    VectorXd res(r.size());
    for(uint32_t i = 0; i < r.size(); i++)
    {
        res(basis.rotl(i, 1)) = r(i);
    }
    return res;
}

template<class Basis> uint64_t countDim1D(Basis&& basis, bool useU1)
{
    const auto N = basis.getN();
    const auto k = basis.getK();
    const int expk = (k == 0) ? 1 : -1;

    std::unordered_set<VectorXi, hash_vector> basisVecs;

    for(uint32_t n = 0; n < (1U << N); ++n)
    {
        if(useU1 && (__builtin_popcountll(n) != N / 2))
        {
            continue;
        }
        Eigen::VectorXi v = Eigen::VectorXi::Zero(1UL << N);
        for(uint32_t k = 0; k < N; k++)
        {
            v(basis.rotl(n, k)) += powi(expk, k);
        }
        if(v.cwiseAbs().sum() != 0)
        {
            make_first_positive(v);
            basisVecs.emplace(std::move(v));
        }
    }

    return basisVecs.size();
}

template<class Basis> uint64_t countDim1DZ2(Basis&& basis, bool useU1)
{
    const auto N = basis.getN();
    const auto k = basis.getK();
    const int p = basis.getP();
    const int expk = (k == 0) ? 1 : -1;

    std::unordered_set<VectorXi, hash_vector> basisVecs;

    for(uint32_t n = 0; n < (1U << N); ++n)
    {
        if(useU1 && (__builtin_popcountll(n) != N / 2))
        {
            continue;
        }
        Eigen::VectorXi v = Eigen::VectorXi::Zero(1UL << N);
        for(uint32_t k = 0; k < N; k++)
        {
            v(basis.rotl(n, k)) += powi(expk, k);
            v(basis.rotl(basis.flip(n), k)) += powi(expk, k) * p;
        }
        if(v.cwiseAbs().sum() != 0)
        {
            make_first_positive(v);
            basisVecs.emplace(std::move(v));
        }
    }

    return basisVecs.size();
}

template<class Basis> void CheckBasis1DParity(Basis&& basis, const MatrixXd& r)
{
    const auto k = basis.getK();
    const int p = basis.getP();

    const int expk = (k == 0) ? 1 : -1;

    for(int i = 0; i < r.cols(); i++)
    {
        VectorXd c = r.col(i);
        REQUIRE((c - expk * translate(basis, c)).norm() < 1e-10);
        REQUIRE((c - p * flip(basis, c)).norm() < 1e-10);
    }
}

template<class Basis> void CheckBasis1D(Basis&& basis, const MatrixXd& r)
{
    const auto k = basis.getK();
    const int expk = (k == 0) ? 1 : -1;

    for(int i = 0; i < r.cols(); i++)
    {
        auto c = r.col(i);
        REQUIRE((c - expk * translate(basis, c)).norm() < 1e-10);
    }
}

template<class Basis> void CheckBasisXXZ(Basis&& basis, const MatrixXd& r)
{
    TIXXZ<uint32_t> tiXXZ(basis, 1.0, 0.9);
    const size_t dim = basis.getDim();
    auto hamTI = edp::constructSparseMat<double>(
        dim, [&tiXXZ](uint32_t col) { return tiXXZ.getCol(col); });

    const auto N = basis.getN();

    XXZ xxz(N, 1.0, 0.9);
    auto hamFull = edp::constructSparseMat<double>(1U << N, xxz);
    Eigen::MatrixXd mat = r.transpose() * hamFull * r;

    REQUIRE((Eigen::MatrixXd(hamTI) - mat).norm() < 1e-10);
}

TEST_CASE("Test Basis1D", "[basis1d]")
{
    SECTION("Representatives are correct")
    {
        uint32_t K = 0;
        for(uint32_t N = 4; N <= 20; N += 4)
        {
            Basis1D<uint32_t> basis(N, K, true);

            for(uint32_t i = 0; i < basis.getDim(); ++i)
            {
                uint32_t rep = basis.getNthRep(i);
                const auto [rep2, rot] = basis.findMinRots(rep);
                REQUIRE(rep2 == rep);
            }
        }
    }
    SECTION("Without U1")
    {
        for(uint32_t N = 4; N <= max_n; N += 2)
        {
            for(uint32_t K : {0U, N / 2})
            {
                DYNAMIC_SECTION("Testing N: " << N << ", K: " << K)
                {
                    Basis1D<uint32_t> basis(N, K, false);

                    REQUIRE(basis.getDim() == countDim1D(basis, false));

                    MatrixXd r = basisMatrix(basis);
                    TestBasisMatrix(r);
                    CheckBasis1D(basis, r);
                    CheckBasisXXZ(basis, r);
                }
            }
        }
    }
    SECTION("With U1")
    {
        for(uint32_t N = 4; N <= max_n; N += 2)
        {
            for(uint32_t K : {0U, N / 2})
            {
                DYNAMIC_SECTION("Testing N: " << N << ", K: " << K)
                {
                    Basis1D<uint32_t> basis(N, K, true);

                    REQUIRE(basis.getDim() == countDim1D(basis, true));

                    MatrixXd r = basisMatrix(basis);
                    TestBasisMatrix(r);
                    CheckBasis1D(basis, r);
                    CheckBasisXXZ(basis, r);
                }
            }
        }
    }
}

TEST_CASE("Test Basis1DZ2", "[basis1dz2]")
{
    const std::array<int, 2> ps{-1, 1};

    SECTION("Without U1")
    {
        for(uint32_t N = 4; N <= max_n; N += 2)
        {
            for(uint32_t K : {0U, N / 2})
            {
                for(int p : ps)
                {
                    DYNAMIC_SECTION("Testing N: " << N << ", K: " << K << ", p: " << p)
                    {
                        Basis1DZ2<uint32_t> basis(N, K, p, false);

                        REQUIRE(basis.getDim() == countDim1DZ2(basis, false));

                        MatrixXd r = basisMatrix(basis);
                        TestBasisMatrix(r);
                        CheckBasis1DParity(basis, r);
                        CheckBasisXXZ(basis, r);
                    }
                }
            }
        }
    }
    SECTION("With U1")
    {
        for(uint32_t N = 6; N <= max_n; N += 2)
        {
            for(uint32_t K : {0U, N / 2})
            {
                for(int p : ps)
                {
                    DYNAMIC_SECTION("Testing N: " << N << ", K: " << K << ", p: " << p)
                    {
                        Basis1DZ2<uint32_t> basis(N, K, p, true);

                        REQUIRE(basis.getDim() == countDim1DZ2(basis, true));

                        MatrixXd r = basisMatrix(basis);
                        TestBasisMatrix(r);
                        CheckBasis1D(basis, r);
                        CheckBasisXXZ(basis, r);
                    }
                }
            }
        }
    }
}
