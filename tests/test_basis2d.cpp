#include "edlib/Basis/Basis2D.hpp"
#include "edlib/Basis/Basis2DZ2.hpp"
#include "edlib/Hamiltonians/TITFI2D.hpp"
#include "edlib/Op/NodeMV.hpp"

#include "edlib/EDP/ConstructSparseMat.hpp"
#include "edlib/EDP/LocalHamiltonian.hpp"

#include "utils.hpp"

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <catch2/catch.hpp>

#include <iostream>
#include <unordered_set>

using namespace edlib;
using namespace Eigen;

template<typename Basis> Eigen::VectorXd translateY(Basis&& basis, const Eigen::VectorXd& r)
{
    Eigen::VectorXd res(r.size());
    for(int i = 0; i < r.size(); i++)
    {
        res(basis.rotateY(i, 1)) = r(i);
    }
    return res;
}

template<typename Basis> Eigen::VectorXd translateX(Basis&& basis, const Eigen::VectorXd& r)
{
    Eigen::VectorXd res(r.size());
    for(int i = 0; i < r.size(); i++)
    {
        res(basis.rotateX(i, 1)) = r(i);
    }
    return res;
}

template<class Basis> void CheckBasis2DParity(Basis&& basis, const MatrixXd& r)
{
    const uint32_t kx = basis.getKx();
    const uint32_t ky = basis.getKy();
    const uint32_t p = basis.getP();

    const int expkx = (kx == 0) ? 1 : -1;
    const int expky = (ky == 0) ? 1 : -1;

    for(int i = 0; i < r.cols(); i++)
    {
        VectorXd c = r.col(i);
        REQUIRE((c - expkx * translateX(basis, c)).norm() < 1e-10);
        REQUIRE((c - expky * translateY(basis, c)).norm() < 1e-10);
        REQUIRE((c - p * flip(basis, c)).norm() < 1e-10);
    }
}

template<class Basis> void CheckBasis2D(Basis&& basis, const MatrixXd& r)
{
    const uint32_t kx = basis.getKx();
    const uint32_t ky = basis.getKy();

    const int expkx = (kx == 0) ? 1 : -1;
    const int expky = (ky == 0) ? 1 : -1;

    for(int i = 0; i < r.cols(); i++)
    {
        VectorXd c = r.col(i);
        REQUIRE((c - expkx * translateX(basis, c)).norm() < 1e-10);
        REQUIRE((c - expky * translateY(basis, c)).norm() < 1e-10);
    }
}

template<class Basis> uint64_t countDim2D(Basis&& basis, bool useU1)
{
    const uint32_t Lx = basis.getLx();
    const uint32_t Ly = basis.getLy();
    const uint32_t kx = basis.getKx();
    const uint32_t ky = basis.getKy();
    const uint32_t N = Lx * Ly;

    const int expkx = (kx == 0) ? 1 : -1;
    const int expky = (ky == 0) ? 1 : -1;

    std::unordered_set<VectorXi, hash_vector> basisVecs;

    for(uint32_t n = 0; n < (1U << N); ++n)
    {
        if(useU1 && (__builtin_popcountll(n) != N / 2))
        {
            continue;
        }
        Eigen::VectorXi v = Eigen::VectorXi::Zero(1UL << N);
        for(uint32_t rx = 0; rx < Lx; rx++)
        {
            for(uint32_t ry = 0; ry < Ly; ry++)
            {
                auto s = basis.rotateX(n, rx);
                s = basis.rotateY(s, ry);
                v(s) += powi(expkx, rx) * powi(expky, ry);
            }
        }
        if(v.cwiseAbs().sum() != 0)
        {
            make_first_positive(v);
            basisVecs.emplace(std::move(v));
        }
    }

    return basisVecs.size();
}

template<class Basis> uint64_t countDim2DZ2(Basis&& basis, bool useU1)
{
    const uint32_t Lx = basis.getLx();
    const uint32_t Ly = basis.getLy();
    const uint32_t kx = basis.getKx();
    const uint32_t ky = basis.getKy();

    const int p = basis.getParity();

    const uint32_t N = Lx * Ly;

    const int expkx = (kx == 0) ? 1 : -1;
    const int expky = (ky == 0) ? 1 : -1;

    std::unordered_set<VectorXi, hash_vector> basisVecs;

    for(uint32_t n = 0; n < (1U << N); ++n)
    {
        if(useU1 && (__builtin_popcountll(n) != N / 2))
        {
            continue;
        }
        Eigen::VectorXi v = Eigen::VectorXi::Zero(1UL << N);
        for(uint32_t rx = 0; rx < Lx; rx++)
        {
            for(uint32_t ry = 0; ry < Ly; ry++)
            {
                auto s = basis.rotateX(n, rx);
                s = basis.rotateY(s, ry);
                v(s) += powi(expkx, rx) * powi(expky, ry);

                s = basis.rotateX(basis.flip(n), rx);
                s = basis.rotateY(s, ry);
                v(s) += p * powi(expkx, rx) * powi(expky, ry);
            }
        }
        if(v.cwiseAbs().sum() != 0)
        {
            make_first_positive(v);
            basisVecs.emplace(std::move(v));
        }
    }

    return basisVecs.size();
}

void testBasis2D(uint32_t Lx, uint32_t Ly, bool useU1)
{
    std::vector<uint32_t> kxs;
    std::vector<uint32_t> kys;

    kxs.push_back(0);
    kys.push_back(0);

    if(Lx % 2 == 0)
    {
        kxs.push_back(Lx / 2);
    }
    if(Ly % 2 == 0)
    {
        kys.push_back(Ly / 2);
    }

    for(const uint32_t kx : kxs)
    {
        for(const uint32_t ky : kys)
        {
            Basis2D<uint32_t> basis(Lx, Ly, kx, ky, useU1);

            DYNAMIC_SECTION("Test Lx: " << Lx << ", Ly: " << Ly << "kx: " << kx << ", ky" << ky)
            {
                INFO("dim: " << basis.getDim());

                REQUIRE(basis.getDim() == countDim2D(basis, useU1));

                MatrixXd r = basisMatrix(basis);
                TestBasisMatrix(r);
                CheckBasis2D(basis, r);
            }
        }
    }
}

void testBasis2DZ2(uint32_t Lx, uint32_t Ly, bool useU1)
{
    std::vector<uint32_t> kxs;
    std::vector<uint32_t> kys;

    kxs.push_back(0);
    kys.push_back(0);

    if(Lx % 2 == 0)
    {
        kxs.push_back(Lx / 2);
    }
    if(Ly % 2 == 0)
    {
        kys.push_back(Ly / 2);
    }

    for(const uint32_t kx : kxs)
    {
        for(const uint32_t ky : kys)
        {
            Basis2DZ2<uint32_t> basis(Lx, Ly, kx, ky, 1, useU1);
            DYNAMIC_SECTION("Test Lx: " << Lx << ", Ly: " << Ly << "kx: " << kx << ", ky" << ky)
            {
                INFO("parity:" << basis.getParity() << ", dim:" << basis.getDim());

                if(basis.getDim() == 0)
                {
                    continue;
                }

                REQUIRE(basis.getDim() == countDim2DZ2(basis, useU1));

                MatrixXd r = basisMatrix(basis);
                TestBasisMatrix(r);
                CheckBasis2D(basis, r);
            }
        }
    }
    for(const uint32_t kx : kxs)
    {
        for(const uint32_t ky : kys)
        {
            Basis2DZ2<uint32_t> basis(Lx, Ly, kx, ky, -1, useU1);
            DYNAMIC_SECTION("Test Lx: " << Lx << ", Ly: " << Ly << "kx: " << kx << ", ky" << ky)
            {
                INFO("parity:" << basis.getParity() << ", dim:" << basis.getDim());

                if(basis.getDim() == 0)
                {
                    continue;
                }

                REQUIRE(basis.getDim() == countDim2DZ2(basis, useU1));

                MatrixXd r = basisMatrix(basis);
                TestBasisMatrix(r);
                CheckBasis2D(basis, r);
            }
        }
    }
}

TEST_CASE("Check Basis2D", "[basis2d]")
{
    SECTION("Not use U1")
    {
        testBasis2D(2, 2, false);
        testBasis2D(3, 2, false);
        testBasis2D(3, 3, false);
        testBasis2D(4, 3, false);
        // testBasis2D(4, 4, false);
        // testBasis2D(5, 4, false);
    }

    SECTION("Use U1")
    {
        testBasis2D(2, 2, true);
        testBasis2D(3, 2, true);
        testBasis2D(3, 3, true);
        testBasis2D(4, 3, true);
        // testBasis2D(4, 4, true);
        // testBasis2D(5, 4, true);
    }
}

TEST_CASE("Check Basis2DZ2Z2", "[basis2dz2]")
{
    SECTION("Not use U1")
    {
        testBasis2DZ2(2, 2, false);
        testBasis2DZ2(3, 2, false);
        testBasis2DZ2(3, 3, false);
        testBasis2DZ2(4, 3, false);
        testBasis2DZ2(4, 4, false);
        // testBasis2DZ2(5, 4, false);
    }

    SECTION("Use U1")
    {
        testBasis2DZ2(2, 2, true);
        testBasis2DZ2(3, 2, true);
        // U(1) symmetry cannot be imposed to 3x3 lattice
        testBasis2DZ2(4, 3, true);
        testBasis2DZ2(4, 4, true);
        // testBasis2DZ2(5, 4, true);
    }
}

TEST_CASE("Compare enegies from the 2D TFI model", "[tfi-2d]")
{
    using Spectra::SortRule;
    using Spectra::CompInfo;
    constexpr uint32_t max_iter = 1000;
    constexpr double eps = 1e-10;

    const uint32_t Lx = 3;
    const uint32_t Ly = 2;
    const uint32_t N = Lx * Ly;
    Basis2D<uint32_t> basis2D(Lx, Ly, 0, 0, false);

    double h = 10.0;

    double gsEnergy1 = 0.0;
    {
        edp::LocalHamiltonian<double> lh(N, 2);
        for(uint32_t nx = 0; nx < Lx; ++nx)
        {
            for(uint32_t ny = 0; ny < Ly; ++ny)
            {
                lh.addTwoSiteTerm(std::make_pair(Lx * ny + nx, Lx * ((ny + 1) % Ly) + nx),
                                  -getSZZ());
                lh.addTwoSiteTerm(std::make_pair(Lx * ny + nx, Lx * ny + ((nx + 1) % Lx)),
                                  -getSZZ());
                lh.addOneSiteTerm(Lx * ny + nx, -h * getSX());
            }
        }
        Eigen::SparseMatrix<double> m = edp::constructSparseMat<double>(1U << N, lh);

        Spectra::SparseSymMatProd<double> op(m);
        Spectra::SymEigsSolver<decltype(op)> eigs(op, 2, 6);
        eigs.init();
        eigs.compute(SortRule::SmallestAlge, max_iter, eps, SortRule::SmallestAlge);
        INFO(eigs.eigenvalues());
        gsEnergy1 = eigs.eigenvalues()[0];
    }

    double gsEnergy2 = 0.0;
    {
        TITFI2D<uint32_t> ham(basis2D, 1.0, h);
        const auto dim = basis2D.getDim();
        NodeMV mv(dim, 0, dim, ham);
        Spectra::SymEigsSolver<NodeMV> eigs(mv, 2, 4);
        eigs.init();
        eigs.compute(SortRule::SmallestAlge, max_iter, eps, SortRule::SmallestAlge);
        INFO(eigs.eigenvalues());
        gsEnergy2 = eigs.eigenvalues()[0];
    }

    REQUIRE_THAT(gsEnergy1, Catch::Matchers::WithinAbs(gsEnergy2, 1e-8));
}
