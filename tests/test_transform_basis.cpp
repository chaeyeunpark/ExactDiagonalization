#include "utils.hpp"

#include "edlib/Basis/Basis1D.hpp"
#include "edlib/Basis/Basis1DZ2.hpp"
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

TEST_CASE("toOriginalBasis", "[TransformBasis]")
{
    std::mt19937 re{1337};
    const uint32_t N = 10;

    SECTION("Test using Basis1D")
    {
        for(const auto& k : {0U, N / 2})
        {
            Basis1D<uint32_t> basis(N, k, false);
            const auto random_vector_before = randomVector(re, basis.getDim());
            const auto random_vector_after
                = toReducedBasis(basis, toOriginalBasis(basis, random_vector_before));

            REQUIRE_THAT(random_vector_before, ::Approx(random_vector_after).epsilon(1e-6));
        }
    }
}
