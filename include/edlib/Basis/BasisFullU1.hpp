#pragma once
#include "AbstractBasis.hpp"
#include "BasisJz.hpp"

#include <tbb/concurrent_vector.h>
#include <tbb/tbb.h>

namespace edlib
{
template<typename UINT> class BasisFullU1 final : public AbstractBasis<UINT>
{
private:
    std::vector<UINT> basis_;

public:
    explicit BasisFullU1(uint32_t N) : AbstractBasis<UINT>(N)
    {
        const uint32_t n = this->getN();
        const uint32_t nup = n / 2;

        BasisJz<UINT> basis(n, nup);

        std::copy(basis.begin(), basis.end(), std::back_inserter(basis_));
    }

    [[nodiscard]] std::size_t getDim() const override { return basis_.size(); }

    [[nodiscard]] UINT getNthRep(uint32_t n) const override { return basis_[n]; }

    [[nodiscard]] auto hamiltonianCoeff(UINT bSigma, [[maybe_unused]] int aidx) const
        -> std::pair<int, double> override
    {
        const auto biter = std::lower_bound(basis_.begin(), basis_.end(), bSigma);

        if(biter == basis_.end() || biter->second != bSigma)
        {
            return {-1, 0.0};
        }

        int bidx = std::distance(basis_.begin(), biter);
        return {bidx, 1.0};
    }

    std::vector<std::pair<UINT, double>> basisVec(uint32_t n) const override
    {
        return {{basis_[n], 1.0}};
    }
};
} // namespace edlib
