#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

#include "AbstractBasis2D.hpp"
#include "BasisJz.hpp"

#include <tbb/tbb.h>

namespace edlib
{
template<typename UINT> class Basis2D final : public AbstractBasis2D<UINT>
{
private:
    // Representatives and number of repretations
    tbb::concurrent_vector<std::pair<UINT, uint32_t>> rpts_;

    void constructBasis(bool useU1)
    {
        if(useU1)
        {
            const uint32_t n = this->getN();
            const uint32_t nup = n / 2;

            BasisJz<UINT> basis(n, nup);

            tbb::parallel_for_each(basis.begin(), basis.end(),
                                   [&](UINT s)
                                   {
                                       uint32_t numRep = this->checkState(s);
                                       if(numRep > 0)
                                       {
                                           rpts_.emplace_back(s, numRep);
                                       }
                                   });
        }
        else
        {
            // if the MSB is 1, the fliped one is smaller.
            auto Lx = this->Lx_;
            auto Ly = this->Ly_;
            tbb::parallel_for(UINT(0U), (UINT(1U) << UINT(Lx * Ly)),
                              [&](UINT s)
                              {
                                  uint32_t numRep = this->checkState(s);
                                  if(numRep > 0)
                                  {
                                      rpts_.emplace_back(s, numRep);
                                  }
                              });
        }

        tbb::parallel_sort(rpts_.begin(), rpts_.end());
    }

public:
    Basis2D(uint32_t Lx, uint32_t Ly, uint32_t kx, uint32_t ky, bool useU1)
        : AbstractBasis2D<UINT>(Lx, Ly, kx, ky)
    {
        assert((!useU1) || (this->getN() % 2 == 0));
        constructBasis(useU1);
    }

    Basis2D(const Basis2D&) = default;
    Basis2D(Basis2D&&) = default;
    Basis2D& operator=(const Basis2D&) = default;
    Basis2D& operator=(Basis2D&&) = default;

    ~Basis2D() = default;

    std::size_t getDim() const override { return rpts_.size(); }

    UINT getNthRep(uint32_t n) const override { return rpts_[n].first; }

    uint32_t repetitions(uint32_t n) const { return rpts_[n].second; }

    inline UINT flip(UINT value) const { return ((this->getUps()) ^ value); }

    std::pair<int, double> hamiltonianCoeff(UINT bSigma, int aidx) const override
    {
        using std::abs;
        using std::pow;
        using std::sqrt;

        uint32_t aNumRep = rpts_[aidx].second;
        double Na = aNumRep;

        UINT bRep;
        uint32_t bRotX, bRotY;
        std::tie(bRep, bRotX, bRotY) = this->getMinRots(bSigma);
        auto iter = std::lower_bound(rpts_.begin(), rpts_.end(), std::make_pair(bRep, 0u));
        if(iter == rpts_.end() || iter->first != bRep)
        {
            return std::make_pair(-1, 0.0);
        }

        int idx = std::distance(rpts_.begin(), iter);
        uint32_t bNumRep = iter->second;
        double Nb = bNumRep;

        return std::make_pair(idx, sqrt(Nb / Na) * this->phase(bRotX, bRotY));
    }

    std::vector<std::pair<UINT, double>> basisVec(uint32_t n) const override
    {
        using std::pow;

        std::map<UINT, double> res;
        UINT rep = rpts_[n].first;
        uint32_t numRep = rpts_[n].second;

        const auto Lx = this->Lx_;
        const auto Ly = this->Ly_;

        double norm = 1.0 / sqrt(numRep * this->getN());

        for(uint32_t nx = 0; nx < Lx; nx++)
        {
            auto srx = this->rotateX(rep, nx);
            for(uint32_t ny = 0; ny < Ly; ny++)
            {
                auto sr = this->rotateY(srx, ny);
                res[sr] += this->phase(nx, ny) * norm;
            }
        }

        return std::vector<std::pair<UINT, double>>(res.begin(), res.end());
    }
};
} // namespace edlib
