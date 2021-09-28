#pragma once

#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cassert>

#include "AbstractBasis2D.hpp"
#include "BasisJz.hpp"

#include <tbb/tbb.h>


namespace edlib
{
template<typename UINT>
class Basis2DZ2
	: public AbstractBasis2D<UINT>
{
public:
	struct RepData
	{
		uint32_t numRep;
		int parity; //Z2 parity
	};

private:
	const int parity_; //Z2_parity

	tbb::concurrent_vector<std::pair<UINT, RepData>> rpts_; //Representatives

	void constructBasis(bool useU1)
	{
		tbb::concurrent_vector<std::tuple<UINT, uint32_t>> candids;
		const unsigned int N = this->getN();
		candids.reserve((1<<(N-3))/N);

		if(useU1)
		{
			const unsigned int n = this->getN();
			const unsigned int nup = n/2;

			BasisJz<UINT> basis(n, nup);

			tbb::parallel_for_each(basis.begin(), basis.end(), [&](UINT s)
			{
				uint32_t numRep = this->checkState(s);
				if(numRep > 0)
				{
					candids.emplace_back(s, numRep);
				}
			});
		}
		else
		{
			//if the MSB is 1, the fliped one is smaller.
			const auto Lx = this->Lx_;
			const auto Ly = this->Ly_;
			tbb::parallel_for(UINT(0u), (UINT(1u) << UINT(Lx*Ly - 1)),
			[&](UINT s) {
				uint32_t numRep = this->checkState(s);
				if(numRep > 0)
				{
					candids.emplace_back(s, numRep);
				}
			});
		}

		tbb::parallel_for(UINT(0u), UINT(candids.size()), 
				[&] (UINT idx) {
			UINT rep;
			uint32_t numRep;
			std::tie(rep, numRep)= candids[idx];
			UINT fliped = flip(rep);

			UINT flipedRep;
			uint32_t rotX, rotY;
			std::tie(flipedRep, rotX, rotY) = this->getMinRots(fliped);

			if(flipedRep == rep && this->phase(rotX, rotY)*parity_ == 1)
			{
				rpts_.emplace_back(rep, RepData{numRep, 0});
			}
			else if(flipedRep > rep)
			{
				rpts_.emplace_back(rep, RepData{numRep, 1});
			}
			else {
				;
			}
		});

		//sort to make it consistent over different instances

		auto comp = [](const std::pair<UINT, RepData>& v1, const std::pair<UINT, RepData>& v2)
		{
			return v1.first < v2.first;
		};

		tbb::parallel_sort(rpts_, comp);
	}


public:
	Basis2DZ2(uint32_t Lx, uint32_t Ly, uint32_t kx, uint32_t ky, 
			int parity, bool useU1)
		: AbstractBasis2D<UINT>(Lx, Ly, kx, ky), parity_(parity)
	{
		assert((parity == 1) || (parity == -1));
		constructBasis(useU1);
	}

	Basis2DZ2<UINT>(const Basis2DZ2<UINT>& ) = default;
	Basis2DZ2<UINT>(Basis2DZ2<UINT>&& ) = default;

	unsigned int stateIdx(UINT rep) const
	{
		auto comp = [](const std::pair<UINT, RepData>& v1, UINT v2)
		{
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

	inline int getParity() const { return parity_; }

	std::size_t getDim() const override
	{
		return rpts_.size();
	}

	UINT getNthRep(int n) const override
	{
		return rpts_[n].first;
	}

	inline UINT flip(UINT value) const
	{
		return ((this->getUps())^value);
	}

	std::pair<int, double> hamiltonianCoeff(UINT bSigma, int aidx) const override
	{
		using std::abs;
		using std::sqrt;
		using std::pow;

		double c = 1.0;
		auto pa = rpts_[aidx].second;
		double Na = pa.numRep/double(1 + pa.parity);

		UINT bRep;
		uint32_t bRotX, bRotY;
		std::tie(bRep, bRotX, bRotY) = this->getMinRots(bSigma);
		auto bidx = stateIdx(bRep);
		if(bidx == getDim())
		{
			c *= parity_;
			std::tie(bRep, bRotX, bRotY) = this->getMinRots(flip(bSigma));
			auto bidx = stateIdx(bRep);

			if(bidx == getDim())
				return std::make_pair(-1, 0.0);
		}

		auto pb = rpts_[bidx].second;
		double Nb = pb.numRep/double(1 + pb.parity);

		return std::make_pair(bidx, c*sqrt(Nb/Na)*this->phase(bRotX,bRotY));
	}

	std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const override
	{
		using std::pow;

		std::map<UINT, double> res;
		UINT rep = getNthRep(n);
		auto repData = rpts_[n].second;
	
		double norm = 1.0/sqrt(repData.numRep*(1 + repData.parity)*this->getN());

		const auto Lx = this->Lx_;
		const auto Ly = this->Ly_;

		for(uint32_t nx = 0; nx < Lx; nx++)
		{
			auto srx = this->rotateX(rep, nx);
			for(uint32_t ny = 0; ny < Ly; ny++)
			{
				auto sr = this->rotateY(srx, ny);
				res[sr] += this->phase(nx, ny)*norm;
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
					res[sr] += this->phase(nx, ny)*parity_*norm;
				}
			}
		}

		return std::vector<std::pair<UINT, double>>(res.begin(), res.end());
	}
};
} // namespace edlib
