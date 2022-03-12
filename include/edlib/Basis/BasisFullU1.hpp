#pragma once
#include "AbstractBasis.hpp"

#include <tbb/tbb.h>

namespace edlib
{
template<typename UINT> class BasisFullU1 final : public AbstractBasis<UINT>
{
private:
public:
    BasisFullU1(uint32_t N) : AbstractBasis<UINT>(N) { }

    std::size_t getDim() const override { return (1u << this->getN()); }

    UINT getNthRep(uint32_t n) const override { return UINT(n); }

    std::pair<int, double> hamiltonianCoeff(UINT bsigma, int aidx) const override
    {
        return std::make_pair(int(bsigma), 1.0);
    }

    std::vector<std::pair<UINT, double>> basisVec(uint32_t n) const override
    {
        return std::vector<std::pair<UINT, double>>{std::make_pair(UINT(n), 1.0)};
    }
};
} // namespace edlib
