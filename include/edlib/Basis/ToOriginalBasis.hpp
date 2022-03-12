#pragma once
#include <tbb/tbb.h>
#include <vector>

namespace edlib
{
/// version with less momory
template<typename UINT, template<typename> class Basis>
std::vector<double> toOriginalVectorLM(const Basis<UINT>& basis, const double* st)
{
    uint32_t N = basis.getN();
    UINT size = UINT(1) << UINT(N);
    std::vector<double> res(size, 0.0);
    tbb::parallel_for(0U, (uint32_t)basis.getDim(), [&](uint32_t n) {
        auto v = basis.basisVec(n);
        for(const auto p : v)
        {
            res[p.first] += st[n] * p.second;
        }
    });
    return res;
}
} // namespace edlib
