#pragma once
#include <mkl.h>
#include <tbb/tbb.h>
#include <vector>

namespace edlib
{
template<typename UINT, template<typename> class Basis> struct SumVector
{
    const Basis<UINT>& basis_;
    const double* st_;
    std::vector<double> my_;

    SumVector(SumVector& x, [[maybe_unused]] tbb::split splt)
        : basis_(x.basis_), st_(x.st_), my_((UINT(1U) << UINT(basis_.getN())), 0.0)
    { }

    SumVector(const Basis<UINT>& basis, const double* st)
        : basis_(basis), st_(st), my_((UINT(1U) << UINT(basis_.getN())), 0.0)
    { }

    void operator()(const tbb::blocked_range<std::uint32_t>& r)
    {
        for(std::uint32_t i = r.begin(); i != r.end(); ++i)
        {
            auto v = basis_.basisVec(i);
            for(const auto p : v)
            {
                my_[p.first] += st_[i] * p.second;
            }
        }
    }
    void join(const SumVector& rhs)
    {
        UINT size = UINT(1) << UINT(basis_.getN());
        double a = 1.0;
        cblas_daxpy(size, a, rhs.my_.data(), 1, my_.data(), 1);
    }
};

template<typename UINT, template<typename> class Basis>
std::vector<double> toOriginalVector(const Basis<UINT>& basis, const double* st)
{
    SumVector<UINT, Basis> sv(basis, st);
    tbb::parallel_reduce(tbb::blocked_range<uint32_t>(0U, basis.getDim()), sv);
    return sv.my_;
}

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