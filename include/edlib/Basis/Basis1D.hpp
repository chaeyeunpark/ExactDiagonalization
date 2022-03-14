#pragma once
#include <cassert>
#include <cmath>
#include <iostream>

#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_sort.h>

#include "AbstractBasis1D.hpp"
#include "BasisJz.hpp"

namespace edlib
{
template<typename UINT> class Basis1D final : public AbstractBasis1D<UINT>
{
private:
    tbb::concurrent_vector<std::pair<UINT, uint32_t>> rpts_; // Representatives

    void constructBasisFull()
    {
        // insert 0
        {
            UINT s = 0;
            int r = this->checkState(s);
            if(r > 0)
            {
                rpts_.emplace_back(s, r);
            }
        }

        // iterate over all odd numbers
        const uint32_t N = this->getN();
        tbb::parallel_for(UINT(1), (UINT(1) << UINT(N)), UINT(2), [&](UINT s) {
            int r = this->checkState(s);
            if(r > 0)
            {
                rpts_.emplace_back(s, r);
            }
        });
        tbb::parallel_sort(rpts_.begin(), rpts_.end());
    }

    void constructBasisJz()
    {
        const uint32_t n = this->getN();
        const uint32_t nup = n / 2;

        BasisJz<UINT> basis(n, nup);

        tbb::parallel_for_each(basis.begin(), basis.end(), [&](UINT s) {
            int r = this->checkState(s);
            if(r > 0)
            {
                rpts_.emplace_back(s, r);
            }
        });
        tbb::parallel_sort(rpts_.begin(), rpts_.end());
    }

public:
    Basis1D(uint32_t N, uint32_t k, bool useU1) : AbstractBasis1D<UINT>(N, k)
    {
        assert((!useU1) || (N % 2 == 0));
        if(useU1)
        {
            constructBasisJz();
        }
        else
        {
            constructBasisFull();
        }
    }

    [[nodiscard]] auto stateIdx(UINT rep) const -> uint32_t
    {
        auto comp = [](const std::pair<UINT, uint32_t>& v1, UINT v2) {
            return v1.first < v2;
        };
        auto iter = lower_bound(rpts_.begin(), rpts_.end(), rep, comp);
        if((iter == rpts_.end()) || (iter->first != rep))
        {
            return getDim();
        }
        else
        {
            return distance(rpts_.begin(), iter);
        }
    }

    [[nodiscard]] auto
    getRepresentatives() const& -> const tbb::concurrent_vector<std::pair<UINT, uint32_t>>&
    {
        return rpts_;
    }
    [[nodiscard]] auto getRepresentatives() && -> tbb::concurrent_vector<std::pair<UINT, uint32_t>>
    {
        return rpts_;
    }

    [[nodiscard]] auto getDim() const -> std::size_t override { return rpts_.size(); }

    [[nodiscard]] auto getNthRep(uint32_t n) const -> UINT override { return rpts_[n].first; }

    [[nodiscard]] inline auto rotRpt(uint32_t n) const -> uint32_t { return rpts_[n].second; }

    [[nodiscard]] auto hamiltonianCoeff(UINT bSigma, int aidx) const
        -> std::pair<int, double> override
    {
        using std::pow;
        using std::sqrt;

        const auto k = this->getK();
        double expk = (k == 0) ? 1.0 : -1.0;

        const auto [bRep, bRot] = this->findMinRots(bSigma);

        auto bidx = stateIdx(bRep);

        if(bidx >= getDim())
        {
            return std::make_pair(-1, 0.0);
        }

        double Na = 1.0 / rpts_[aidx].second;
        double Nb = 1.0 / rpts_[bidx].second;

        return std::make_pair(bidx, sqrt(Nb / Na) * pow(expk, bRot));
    }

    [[nodiscard]] auto basisVec(uint32_t n) const -> std::vector<std::pair<UINT, double>> override
    {
        const auto k = this->getK();
        const double expk = (k == 0) ? 1.0 : -1.0;
        std::vector<std::pair<UINT, double>> res;

        auto rep = rpts_[n].first;
        double norm = 1.0 / sqrt(rpts_[n].second);
        for(uint32_t r = 0; r < rpts_[n].second; r++)
        {
            res.emplace_back(this->rotl(rep, r), pow(expk, r) * norm);
        }
        return res;
    }
};
} // namespace edlib
