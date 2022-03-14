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
template<typename UINT> class Basis2DZ2 final : public AbstractBasis2D<UINT>
{
public:
    struct RepData
    {
        uint32_t numRep;
        int parity; // Z2 parity
    };

private:
    const int parity_; // Z2_parity

    tbb::concurrent_vector<std::pair<UINT, RepData>> rpts_; // Representatives

    void constructBasis(bool useU1)
    {
        tbb::concurrent_vector<std::tuple<UINT, uint32_t>> candids;
        const uint32_t N = this->getN();
        candids.reserve((1U << (N - 3)) / N);

        if(useU1)
        {
            const uint32_t n = this->getN();
            const uint32_t nup = n / 2;

            BasisJz<UINT> basis(n, nup);

            tbb::parallel_for_each(basis.begin(), basis.end(), [&](UINT s) {
                uint32_t numRep = this->checkState(s);
                if(numRep > 0)
                {
                    candids.emplace_back(s, numRep);
                }
            });
        }
        else
        {
            // if the MSB is 1, the fliped one is smaller.
            const auto Lx = this->Lx_;
            const auto Ly = this->Ly_;
            tbb::parallel_for(UINT(0U), (UINT(1U) << UINT(Lx * Ly - 1)), [&](UINT s) {
                uint32_t numRep = this->checkState(s);
                if(numRep > 0)
                {
                    candids.emplace_back(s, numRep);
                }
            });
        }

        tbb::parallel_for(UINT(0U), UINT(candids.size()), [&](UINT idx) {
            const auto [rep, numRep] = candids[idx];
            UINT fliped = flip(rep);

            const auto [flipedRep, rotX, rotY] = this->getMinRots(fliped);

            if(flipedRep == rep && this->phase(rotX, rotY) * parity_ == 1)
            {
                rpts_.emplace_back(rep, RepData{numRep, 0});
            }
            else if(flipedRep > rep)
            {
                rpts_.emplace_back(rep, RepData{numRep, 1});
            }
            else
            {
                ;
            }
        });

        // sort to make it consistent over different instances
        const auto comp = [](const std::pair<UINT, RepData>& v1, const std::pair<UINT, RepData>& v2) {
            return v1.first < v2.first;
        };

        tbb::parallel_sort(rpts_, comp);
    }

public:
    Basis2DZ2(uint32_t Lx, uint32_t Ly, uint32_t kx, uint32_t ky, int parity, bool useU1)
        : AbstractBasis2D<UINT>(Lx, Ly, kx, ky), parity_(parity)
    {
        assert((parity == 1) || (parity == -1));
        constructBasis(useU1);
    }

    [[nodiscard]] auto stateIdx(UINT rep) const -> uint32_t
    {
        const auto comp = [](const std::pair<UINT, RepData>& v1, UINT v2) {
            return v1.first < v2;
        };
        const auto iter = lower_bound(rpts_.begin(), rpts_.end(), rep, comp);
        if((iter == rpts_.end()) || (iter->first != rep))
        {
            return getDim();
        }
        else
        {
            return distance(rpts_.begin(), iter);
        }
    }

    [[nodiscard]] inline auto getParity() const -> int { return parity_; }

    [[nodiscard]] auto getDim() const -> std::size_t override { return rpts_.size(); }

    [[nodiscard]] UINT getNthRep(uint32_t n) const override { return rpts_[n].first; }

    [[nodiscard]] inline UINT flip(UINT value) const { return ((this->getUps()) ^ value); }

    [[nodiscard]] auto hamiltonianCoeff(UINT bSigma, int aidx) const
        -> std::pair<int, double> override
    {
        using std::abs;
        using std::pow;
        using std::sqrt;

        double c = 1.0;
        auto pa = rpts_[aidx].second;
        double Na = pa.numRep / double(1 + pa.parity);

        auto [bRep, bRotX, bRotY] = this->getMinRots(bSigma);
        auto bidx = stateIdx(bRep);
        if(bidx == getDim())
        {
            c *= parity_;
            std::tie(bRep, bRotX, bRotY) = this->getMinRots(flip(bSigma));
            bidx = stateIdx(bRep);

            if(bidx == getDim())
            {
                return std::make_pair(-1, 0.0);
            }
        }

        auto pb = rpts_[bidx].second;
        double Nb = pb.numRep / double(1 + pb.parity);

        return std::make_pair(bidx, c * sqrt(Nb / Na) * this->phase(bRotX, bRotY));
    }

    [[nodiscard]] auto basisVec(uint32_t n) const
        -> std::vector<std::pair<UINT, double>> override
    {
        using std::pow;

        std::map<UINT, double> res;
        UINT rep = getNthRep(n);
        auto repData = rpts_[n].second;

        double norm = 1.0 / sqrt(repData.numRep * (1 + repData.parity) * this->getN());

        const auto Lx = this->Lx_;
        const auto Ly = this->Ly_;

        for(uint32_t nx = 0; nx < Lx; nx++)
        {
            auto srx = this->rotateX(rep, nx);
            for(uint32_t ny = 0; ny < Ly; ny++)
            {
                auto sr = this->rotateY(srx, ny);
                res[sr] += this->phase(nx, ny) * norm;
            }
        }
        if(repData.parity != 0)
        {
            UINT fliped = flip(rep);
            for(uint32_t nx = 0; nx < Lx; nx++)
            {
                auto srx = this->rotateX(fliped, nx);
                for(uint32_t ny = 0; ny < Ly; ny++)
                {
                    auto sr = this->rotateY(srx, ny);
                    res[sr] += this->phase(nx, ny) * parity_ * norm;
                }
            }
        }

        return {res.begin(), res.end()};
    }
};
} // namespace edlib
