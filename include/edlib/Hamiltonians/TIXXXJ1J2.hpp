#pragma once
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <cstdint>
#include <map>
//#include "BitOperations.h"
#include "../Basis/AbstractBasis1D.hpp"

template<typename UINT> class TIXXXJ1J2
{
private:
    const edlib::AbstractBasis<UINT>& basis_;
    double J1_;
    double J2_;

    int sign_ = 1;

public:
    TIXXXJ1J2(const edlib::AbstractBasis<UINT>& basis, double J1, double J2, bool signRule = false)
        : basis_(basis), J1_(J1), J2_(J2)
    {
        if(signRule)
            sign_ = -1;
    }

    std::map<int, double> getCol(UINT n) const
    {
        int N = basis_.getN();

        UINT a = basis_.getNthRep(n);
        const boost::dynamic_bitset<> bs(N, a);

        std::map<int, double> m;
        for(unsigned int i = 0; i < N; i++)
        {
            // Nearest neighbors
            {
                unsigned int j = (i + 1) % N;
                int sgn = (1 - 2 * bs[i]) * (1 - 2 * bs[j]);

                m[n] += J1_ * sgn;

                UINT s = a;
                s ^= basis_.mask({i, j});

                int bidx;
                double coeff;

                std::tie(bidx, coeff) = basis_.hamiltonianCoeff(s, n);

                if(bidx >= 0)
                    m[bidx] += J1_ * (1 - sgn) * sign_ * coeff;
            }
            // Next-nearest neighbors
            {
                unsigned int j = (i + 2) % N;
                int sgn = (1 - 2 * bs[i]) * (1 - 2 * bs[j]);

                m[n] += J2_ * sgn;

                UINT s = a;
                s ^= basis_.mask({i, j});

                int bidx;
                double coeff;

                std::tie(bidx, coeff) = basis_.hamiltonianCoeff(s, n);

                if(bidx >= 0)
                    m[bidx] += J2_ * (1 - sgn) * coeff;
            }
        }
        return m;
    }
};
