#pragma once
#include "AbstractBasis.hpp"
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace edlib
{
template<typename UINT> class AbstractBasis2D : public AbstractBasis<UINT>
{
protected:
    const uint32_t Lx_;
    const uint32_t Ly_;
    const uint32_t kx_;
    const uint32_t ky_;
    const UINT px_ = (~UINT(0)) >> (sizeof(UINT) * CHAR_BIT - Lx_);

    /**
     * @brief Get minimum values among all possible x and y rotations.
     *
     * @param sigma The basis state
     * @return Tuple of the minimum rotation, the number of x rotations,
     * and the number of y rotations.
     * */
    [[nodiscard]] auto getMinRots(UINT sigma) const -> std::tuple<UINT, uint32_t, uint32_t>
    {
        /* In two dimension, there are several different rotations that gives
         * the same representation */
        UINT rep = sigma;
        uint32_t rotX = 0;
        uint32_t rotY = 0;
        for(uint32_t rx = 0; rx < Lx_; rx++)
        {
            UINT srx = rotateX(sigma, rx);
            for(uint32_t ry = 0; ry < Ly_; ry++)
            {
                UINT sr = rotateY(srx, ry);
                if(sr < rep)
                {
                    rep = sr;
                    rotX = rx;
                    rotY = ry;
                }
            }
        }
        return std::make_tuple(rep, rotX, rotY);
    }

    /**
     * @brief Check the given basis state can be a representative.
     *
     * @return How many times \f$s\f$ appears among all rotations.
     */
    [[nodiscard]] uint32_t checkState(UINT s) const
    {
        uint32_t cnt = 0;
        UINT sr = s;

        for(uint32_t rx = 1; rx <= Lx_; rx++)
        {
            for(uint32_t ry = 1; ry <= Ly_; ry++)
            {
                sr = rotateX(s, rx);
                sr = rotateY(sr, ry);
                if(sr < s)
                {
                    return 0; // s is not a representative
                }
                else if(sr == s)
                {
                    if(phase(rx, ry) == -1) // not allowed
                    {
                        return 0;
                    }
                    else
                    {
                        ++cnt;
                    }
                }
            }
        }
        return cnt;
    }

    [[nodiscard]] int phase(int32_t rotX, int32_t rotY) const
    {
        // return exp(2\Pi I(kx_*rotX/Lx_ + ky_*rotY/Ly_))
        int sgn = 1;
        if((kx_ * rotX) % Lx_ != 0)
        {
            sgn *= -1;
        }
        if((ky_ * rotY) % Ly_ != 0)
        {
            sgn *= -1;
        }
        return sgn;
    }

public:
    AbstractBasis2D(uint32_t Lx, uint32_t Ly, uint32_t kx, uint32_t ky)
        : AbstractBasis<UINT>(Lx * Ly), Lx_{Lx}, Ly_{Ly}, kx_{kx}, ky_{ky}
    { }

    [[nodiscard]] UINT rotateY(UINT sigma, uint32_t r) const
    {
        const auto N = this->getN();
        const uint32_t count = Lx_ * r;
        return ((sigma << count) & this->getUps()) | (sigma >> (N - count));
    }

    [[nodiscard]] UINT rotateX(UINT sigma, uint32_t r) const
    {
        assert(r <= Lx_);
        const auto Lx = Lx_;
        const auto Ly = Ly_;
        const auto px = px_;
        for(uint32_t ny = 0; ny < Ly; ny++)
        {
            UINT t = (sigma >> (ny * Lx)) & px;
            t = ((t << r) | (t >> (Lx - r))) & px;
            sigma &= ~(px << ny * Lx);
            sigma |= (t << ny * Lx);
        }
        return sigma;
    }

    [[nodiscard]] inline uint32_t getLx() const { return Lx_; }
    [[nodiscard]] inline uint32_t getLy() const { return Ly_; }
    [[nodiscard]] inline uint32_t getKx() const { return kx_; }
    [[nodiscard]] inline uint32_t getKy() const { return ky_; }

    [[nodiscard]] inline uint32_t toIdx(uint32_t nx, uint32_t ny) const { return ny * Lx_ + nx; }
};
} // namespace edlib
