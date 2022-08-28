#pragma once
/**
 * @file
 * Transform basis
 */
#include <tbb/tbb.h>

#include <cassert>
#include <span>
#include <vector>

namespace edlib
{
/**
 * @brief Convert reduced basis state `st` to a state in the original basis
 *
 * @tparam UINT unsigned int type
 * @tparam Basis basis type
 * @param basis reduced basis
 * @param st State in a reduced basis
 */
template<typename UINT, template<typename> class Basis>
std::vector<double> toOriginalBasis(const Basis<UINT>& basis, const std::span<const double> st)
{
    assert(st.size() == basis.getDim());

    const uint32_t N = basis.getN();
    const UINT size = UINT(1) << UINT(N);
    std::vector<double> res(size, 0.0);
    tbb::parallel_for(0U, (uint32_t)basis.getDim(), [&](uint32_t n) {
        const auto v = basis.basisVec(n);
        for(const auto& p : v)
        {
            res[p.first] += st[n] * p.second;
        }
    });
    return res;
}

/**
 * @brief Convert original basis state to the reduced state
 *
 * This function only outputs correct results when the given state is already within the reduced
 * basis.
 *
 * @tparam UINT unsigned int type
 * @tparam Basis basis type
 * @param basis reduced basis
 * @param st State in a reduced basis
 */
template<typename UINT, template<typename> class Basis>
std::vector<double> toReducedBasis(const Basis<UINT>& basis, const std::span<const double> st)
{
    assert(st.size() == (size_t{1} << size_t{basis.getN()}));
    std::vector<double> res(basis.getDim(), 0.0);
    tbb::parallel_for(0U, (uint32_t)basis.getDim(), [&](uint32_t n) {
        const auto v = basis.basisVec(n);

        double r = 0.0;

        for(const auto& p : v)
        {
            r += st[p.first] * p.second;
        }
        res[n] = r;
    });
    return res;
}

} // namespace edlib
