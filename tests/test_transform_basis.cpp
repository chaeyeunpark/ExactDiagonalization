#include "utils.hpp"

#include "edlib/Basis/Basis1D.hpp"
#include "edlib/Basis/Basis2D.hpp"
#include "edlib/Basis/TransformBasis.hpp"
#include "edlib/EDP/ConstructSparseMat.hpp"

#include <Eigen/Dense>
#include <catch2/catch_all.hpp>

#include <algorithm>
#include <iostream>
#include <random>
#include <unordered_set>

using namespace Eigen;
using namespace edlib;

template<typename RandomEngine> std::vector<double> randomVector(RandomEngine& re, size_t size)
{
    std::normal_distribution<double> ndist;
    std::vector<double> res;
    res.reserve(size);

    for(size_t i = 0; i < size; i++)
    {
        res.push_back(ndist(re));
    }
    const double norm = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
    std::for_each(res.begin(), res.end(),
                  [sqrt_norm = std::sqrt(norm)](double& elt) { elt /= sqrt_norm; });
    return res;
}

template<typename T, typename UINT, typename Derived>
Eigen::Matrix<T, Eigen::Dynamic, 1> toTranslationalInvariance(const Basis1D<UINT>& basis,
                                                              const Eigen::DenseBase<Derived>& vec)
{
    const auto N = basis.getN();
    assert(static_cast<size_t>(vec.size()) == (size_t{1U} << N));

    Eigen::Matrix<T, Eigen::Dynamic, 1> res(1U << N);

    for(uint32_t count = 0; count < N; count++)
    {
        for(uint32_t i = 0; i < (1U << N); i++)
        {
            res.coeffRef(i) += vec.coeff(basis.rotl(i, count));
        }
    }
    res.normalize();
    return res;
}

TEST_CASE("toOriginalBasis", "[TransformBasis]")
{
    std::mt19937 re{1337};

    SECTION("Test using Basis1D")
    {
        const uint32_t N = 10;
        for(const auto& k : {0U, N / 2})
        {
            Basis1D<uint32_t> basis(N, k, false);
            const auto random_vector_before = randomVector(re, basis.getDim());
            const auto random_vector_after
                = toReducedBasis(basis, toOriginalBasis(basis, random_vector_before));

            REQUIRE_THAT(random_vector_before, ::Approx(random_vector_after).epsilon(1e-6));
        }
    }

    SECTION("Test using Basis2D")
    {
        const uint32_t Lx = 4;
        const uint32_t Ly = 3;
        for(const auto& k1 : {0U, Lx / 2})
        {
            Basis2D<uint32_t> basis(Lx, Ly, k1, 0, false);
            const auto random_vector_before = randomVector(re, basis.getDim());
            const auto random_vector_after
                = toReducedBasis(basis, toOriginalBasis(basis, random_vector_before));

            REQUIRE_THAT(random_vector_before, ::Approx(random_vector_after).epsilon(1e-6));
        }
    }
}

TEST_CASE("toReducedBasis", "[TransformBasis]")
{
    const uint32_t N = 10;
    std::mt19937 re{1337};

    for(const auto& k : {0U, N / 2})
    {
        Basis1D<uint32_t> basis(N, k, false);
        const Eigen::MatrixXd isometry = basisMatrix(basis);
        const auto vec = randomVector(re, 1U << N);
        const auto full_vector = toTranslationalInvariance<double>(
            basis, Eigen::Map<const VectorXd>(vec.data(), static_cast<Eigen::Index>(vec.size())));
        const Eigen::VectorXd reduced_vector1 = [&basis, &full_vector]() -> Eigen::VectorXd {
            const auto v = toReducedBasis(
                basis, {full_vector.data(), static_cast<size_t>(full_vector.size())});
            return Eigen::Map<const Eigen::VectorXd>(v.data(), static_cast<Eigen::Index>(v.size()));
        }();
        const Eigen::VectorXd reduced_vector2 = full_vector.transpose() * isometry;

        REQUIRE((reduced_vector1 - reduced_vector2).squaredNorm() < 1e-6);
    }
}
