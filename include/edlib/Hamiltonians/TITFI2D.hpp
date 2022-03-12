#pragma once
#include "dynamic_bitset.hpp"
#include "edlib/Basis/AbstractBasis2D.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <map>

template<typename UINT> class TITFI2D
{
private:
    edlib::AbstractBasis2D<UINT>& basis_;
    double J_;
    double h_;

public:
    TITFI2D(edlib::AbstractBasis2D<UINT>& basis, double J, double h) : basis_(basis), J_(J), h_(h)
    { }

    std::map<std::size_t, double> getCol(UINT n) const
    {
        uint32_t Lx = basis_.getLx();
        uint32_t Ly = basis_.getLy();

        UINT a = basis_.getNthRep(n);

        std::map<std::size_t, double> m;

        const edlib::dynamic_bitset bs(Lx * Ly, a);

        for(uint32_t nx = 0; nx < Lx; ++nx)
        {
            for(uint32_t ny = 0; ny < Ly; ++ny)
            {
                // Next-nearest
                const auto idx = basis_.toIdx(nx, ny);
                const int sgn1 = (1 - 2 * static_cast<int>(bs[idx]))
                                 * (1 - 2 * static_cast<int>(bs[basis_.toIdx((nx + 1) % Lx, ny)]));
                const int sgn2 = (1 - 2 * static_cast<int>(bs[idx]))
                                 * (1 - 2 * static_cast<int>(bs[basis_.toIdx(nx, (ny + 1) % Ly)]));

                m[n] += -J_ * (sgn1 + sgn2); // ZZ

                UINT s = a;
                UINT t = (UINT(1U) << UINT(idx)); // X
                s ^= t;

                const auto [bidx, coeff] = basis_.hamiltonianCoeff(s, n);

                if(bidx >= 0)
                {
                    m[bidx] += -h_ * coeff;
                }
            }
        }
        return m;
    }
};
