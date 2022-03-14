#pragma once
#include "AbstractBasis.hpp"
#include <cassert>
#include <cmath>

namespace edlib
{
template<typename UINT> class BasisFullZ2 final : public AbstractBasis<UINT>
{
private:
    int parity_;

public:
    BasisFullZ2(uint32_t N, int parity) : AbstractBasis<UINT>(N), parity_{parity}
    {
        assert(parity == -1 || parity == 1);
    }

    [[nodiscard]] auto getDim() const -> std::size_t override { return std::size_t(1) << (this->getN() - 1); }

    [[nodiscard]] auto getNthRep(uint32_t n) const -> UINT override { return n; }

    [[nodiscard]]
    auto hamiltonianCoeff(UINT bsigma, [[maybe_unused]] int aidx) const
        -> std::pair<int, double> override
    {
        if(bsigma < flip(bsigma))
        {
            return std::make_pair(int(bsigma), 1.0);
        }
        else
        {
            return std::make_pair(int(flip(bsigma)), parity_);
        }
    }

    [[nodiscard]] auto basisVec(uint32_t n) const
        -> std::vector<std::pair<UINT, double>> override
    {
        using std::sqrt;
        return {
            std::make_pair(UINT(n), 1.0 / sqrt(2.0)),
            std::make_pair(UINT(flip(n)), parity_ * 1.0 / sqrt(2.0)),
        };
    }

    [[nodiscard]] inline auto flip(UINT value) const -> UINT { return ((this->getUps()) ^ value); }
};
} // namespace edlib
