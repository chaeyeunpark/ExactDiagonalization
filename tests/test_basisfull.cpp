#include <iostream>
#include <unordered_set>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <boost/dynamic_bitset.hpp>

#include "edlib/Basis/BasisFull.hpp"
#include "edlib/Basis/BasisFullZ2.hpp"
#include "edlib/Basis/ToOriginalBasis.hpp"
#include "edlib/EDP/ConstructSparseMat.hpp"
#include "edlib/Op/NodeMV.hpp"

#include "XXZ.hpp"
#include "utils.hpp"

#include <catch2/catch.hpp>

using namespace Eigen;
using namespace edlib;

constexpr uint32_t MAX_N = 12;

template<typename UINT> class OpenTFI
{
private:
    const edlib::AbstractBasis<UINT>& basis_;
    double J_;
    double h_;

public:
    OpenTFI(const edlib::AbstractBasis<UINT>& basis, double J, double h)
        : basis_{basis}, J_{J}, h_{h}
    { }

    [[nodiscard]] std::map<int, double> getCol(UINT n) const
    {
        std::map<int, double> m;
        const uint32_t N = basis_.getN();

        UINT a = basis_.getNthRep(n);
        const boost::dynamic_bitset<> bs(N, a);

        for(uint32_t i = 0; i < N - 1; ++i)
        {
            int zz = (1 - 2 * static_cast<int>(bs[i])) * (1 - 2 * static_cast<int>(bs[i + 1]));

            m[n] += -J_ * zz;

            UINT s = a;
            s ^= basis_.mask({i});

            const auto [bidx, coeff] = basis_.hamiltonianCoeff(s, n);

            if(bidx >= 0)
            {
                m[bidx] += -h_ * coeff;
            }
        }
        return m;
    }

    [[nodiscard]] std::map<int, double> operator()(UINT n) const { return getCol(n); }
};

TEST_CASE("Test BasisFull and BasisFullZ2 using the transverse field Ising model", "[basisfull]")
{
    constexpr uint32_t max_iter = 10000;
    constexpr double eps = 1e-12;

    constexpr double h = 0.5;
    for(uint32_t N = 4; N <= MAX_N; N += 2)
    {
        BasisFull<uint32_t> basisFull(N);
        BasisFullZ2<uint32_t> basisFullP(N, 1);
        BasisFullZ2<uint32_t> basisFullM(N, -1);

        auto hamFull = OpenTFI(basisFull, 1.0, h);
        auto hamFullP = OpenTFI(basisFullP, 1.0, h);
        auto hamFullM = OpenTFI(basisFullM, 1.0, h);

        VectorXd evFull;
        VectorXd evFullP;
        VectorXd evFullM;

        {
            const auto dim = basisFull.getDim();

            NodeMV mv(dim, 0, dim, hamFull);

            Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, NodeMV> eigs(&mv, 2, 6);
            eigs.init();
            eigs.compute(max_iter, eps, Spectra::SMALLEST_ALGE);
            if(eigs.info() != Spectra::SUCCESSFUL) {
                REQUIRE(false);
            }
            evFull = eigs.eigenvalues();
        }

        {
            const auto dim = basisFullP.getDim();

            NodeMV mv(dim, 0, dim, hamFullP);

            Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, NodeMV> eigs(&mv, 2, 6);
            eigs.init();
            eigs.compute(max_iter, eps, Spectra::SMALLEST_ALGE);
            if(eigs.info() != Spectra::SUCCESSFUL) {
                REQUIRE(false);
            }
            evFullP = eigs.eigenvalues();
        }

        {
            const auto dim = basisFullM.getDim();

            NodeMV mv(dim, 0, dim, hamFullM);

            Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, NodeMV> eigs(&mv, 2, 6);
            eigs.init();
            eigs.compute(max_iter, eps, Spectra::SMALLEST_ALGE);
            if(eigs.info() != Spectra::SUCCESSFUL) {
                REQUIRE(false);
            }
            evFullM = eigs.eigenvalues();
        }
        using Catch::WithinAbs;
        REQUIRE_THAT(evFull(0), WithinAbs(evFullP(0), 1e-8));
        REQUIRE_THAT(evFull(1), WithinAbs(evFullM(0), 1e-8));
    }
}
