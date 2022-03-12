#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <map>
//#include "BitOperations.h"
#include "../Basis/AbstractBasis1D.hpp"

template<typename UINT> class TITFIsing
{
private:
    const edlib::AbstractBasis1D<UINT>& basis_;
    double J_;
    double h_;

public:
    TITFIsing(const edlib::AbstractBasis1D<UINT>& basis, double J, double h)
        : basis_(basis), J_(J), h_(h)
    { }

    std::map<std::size_t, double> getCol(UINT n) const
    {
        unsigned int N = basis_.getN();

        UINT a = basis_.getNthRep(n);

        std::map<std::size_t, double> m;
        for(unsigned int i = 0; i < N; i++)
        {
            // Next-nearest
            {
                const unsigned int j = (i + 1) % N;
                const int bi = static_cast<int>((a >> i) & 1U);
                const int bj = static_cast<int>((a >> j) & 1U);
                const int sgn = (1 - 2 * bi) * (1 - 2 * bj);

                m[n] += -J_ * sgn; // ZZ

                UINT s = a;
                s ^= basis_.mask({i});

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
