#pragma once
#include "AbstractBasis.hpp"

namespace edlib
{
/**
 * @brief Base class for 1D Hamiltonians.
 */
template<typename UINT> class AbstractBasis1D : public AbstractBasis<UINT>
{
protected:
    const uint32_t k_;

    /**
     * @brief Left rotation required to make a basis state the smallest.
     *
     * @param s The basis state as a number in the binary representation
     *
     * TODO: may change to std::optional
     */
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

    /**
     * @brief Rotate left the value in count times.
     *
     * @param value Value to rotate
     * @param count Number of translations
     */
    UINT rotl(UINT value, uint32_t count) const
    {
        const auto N = this->getN();
        count %= N;
        return ((value << count) & this->getUps()) | (value >> (N - count));
    }

    /**
     * @brief Find the smallest value among all rotations.
     *
     * @param sigma Value of the basis vector in the binary representation
     * @return Smallest value and the corresponding rotation pair
     */
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

    /**
     * @brief Eigenvalue of the translational invariant
     * \f$T| \sigma(k) \rangle = e^{ik}| \sigma(k) \ranlge \f$
     */
    [[nodiscard]] inline uint32_t getK() const { return k_; }
};
} // namespace edlib
