#pragma once
#include "dynamic_bitset.hpp"
#include "edlib/Basis/AbstractBasis2D.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <map>

template<typename UINT> class TIJ1J2_2D
{
private:
    const edlib::AbstractBasis2D<UINT>& basis_;
    double J1_;
    double J2_;
    int sign_;

public:
    TIJ1J2_2D(const edlib::AbstractBasis2D<UINT>& basis, double J1, double J2,
              bool signRule = false)
        : basis_(basis), J1_(J1), J2_(J2)
    {
        if(signRule)
        {
            sign_ = -1;
        }
        else
        {
            sign_ = 1;
        }
    }

    std::map<int, double> getCol(UINT n) const
    {
        const uint32_t Lx = basis_.getLx();
        const uint32_t Ly = basis_.getLy();

        const UINT a = basis_.getNthRep(n);
        const edlib::dynamic_bitset bs(Lx * Ly, a);

        std::map<int, double> m;
        for(uint32_t nx = 0; nx < Lx; nx++)
        {
            for(uint32_t ny = 0; ny < Ly; ny++)
            {
                auto i = basis_.toIdx(nx, ny);

                const auto neighbors_indices
                    = {basis_.toIdx((nx + 1) % Lx, ny), basis_.toIdx(nx, (ny + 1) % Ly),
                       basis_.toIdx((nx + 1) % Lx, (ny + 1) % Ly),
                       basis_.toIdx((nx + 1) % Lx, (ny - 1 + Ly) % Ly)};
                for(const auto j : neighbors_indices)
                {
                    const int zz = (1 - 2 * bs[i]) * (1 - 2 * bs[j]);

                    m[n] += J1_ * zz;

                    UINT s = a;
                    s ^= basis_.mask({i, j});

                    const auto [bidx, coeff] = basis_.hamiltonianCoeff(s, n);

                    if(bidx >= 0)
                    {
                        m[bidx] += J1_ * (1 - zz) * coeff * sign_;
                    }
                }
            }
        }
        return m;
    }
};
