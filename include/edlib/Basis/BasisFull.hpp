#pragma once
#include "AbstractBasis.hpp"

#include <tbb/tbb.h>

namespace edlib
{
template<typename UINT> class BasisFull final : public AbstractBasis<UINT>
{
private:
public:
    explicit BasisFull(uint32_t N) : AbstractBasis<UINT>(N) { }

    [[nodiscard]] auto getDim() const -> std::size_t override { return (1U << this->getN()); }

    [[nodiscard]] auto getNthRep(uint32_t n) const-> UINT override { return UINT(n); }

    [[nodiscard]] auto hamiltonianCoeff(UINT bsigma, [[maybe_unused]] int aidx) const
        -> std::pair<int, double> override
    {
        return std::make_pair(int(bsigma), 1.0);
    }

    [[nodiscard]] auto basisVec(uint32_t n) const -> std::vector<std::pair<UINT, double>> override
    {
        return std::vector<std::pair<UINT, double>>{std::make_pair(UINT(n), 1.0)};
    }
};
} // namespace edlib
