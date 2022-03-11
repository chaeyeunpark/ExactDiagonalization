#pragma once
#include "AbstractBasis.hpp"

namespace edlib
{
template<typename UINT> class AbstractBasis1D : public AbstractBasis<UINT>
{
protected:
    const uint32_t k_;

    // TODO: may change to std::optional
    int checkState(UINT s) const
    {
        UINT sr = s;
        const auto N = this->getN();
        for(uint32_t r = 1; r <= N; r++)
        {
            sr = this->rotl(s, r);
            if(sr < s)
            {
                return -1; // s is not a representative
            }
            else if(sr == s)
            {
                if((k_ % (N / r)) != 0)
                {
                    return -1; // this representative is not allowed for this k
                }
                return static_cast<int>(r);
            }
        }
        return -1;
    }

public:
    AbstractBasis1D(uint32_t N, uint32_t k) : AbstractBasis<UINT>(N), k_{k} { }

    UINT rotl(UINT value, uint32_t count) const
    {
        const auto N = this->getN();
        count %= N;
        return ((value << count) & this->getUps()) | (value >> (N - count));
    }

    std::pair<UINT, uint32_t> findMinRots(UINT sigma) const
    {
        const uint32_t N = this->getN();
        UINT rep = sigma;
        uint32_t rot = 0U;
        for(uint32_t r = 1; r < N; r++)
        {
            UINT sr = rotl(sigma, r);
            if(sr < rep)
            {
                rep = sr;
                rot = r;
            }
        }
        return std::make_pair(rep, rot);
    }

    [[nodiscard]] inline uint32_t getK() const { return k_; }
};
} // namespace edlib
