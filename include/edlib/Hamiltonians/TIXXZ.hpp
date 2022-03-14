#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <map>

#include "../Basis/AbstractBasis1D.hpp"

template<typename UINT> class TIXXZ
{
private:
    const edlib::AbstractBasis1D<UINT>& basis_;
    double J_;
    double delta_;

public:
    TIXXZ(const edlib::AbstractBasis1D<UINT>& basis, double J, double delta)
        : basis_(basis), J_(J), delta_(delta)
    { }

    [[nodiscard]] std::map<int, double> getCol(UINT n) const
    {
        unsigned int N = basis_.getN();

        UINT a = basis_.getNthRep(n);

        std::map<int, double> m;
        for(unsigned int i = 0; i < N; i++)
        {
            // Next-nearest
            {
                const unsigned int j = (i + 1) % N;
                const int bi = static_cast<int>((a >> i) & 1U);
                const int bj = static_cast<int>((a >> j) & 1U);
                const int sgn = (1 - 2 * bi) * (1 - 2 * bj);

                m[n] += J_ * delta_ * sgn;

                UINT s = a;
                s ^= basis_.mask({i, j});

                const auto [bidx, coeff] = basis_.hamiltonianCoeff(s, n);

                if(bidx >= 0)
                {
                    m[bidx] += J_ * (1.0 - sgn) * coeff;
                }
            }
        }
        return m;
    }
};
