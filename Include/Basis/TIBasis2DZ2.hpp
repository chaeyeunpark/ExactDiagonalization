#pragma once

#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cassert>

#include "Basis2D.hpp"
#include "BasisJz.hpp"

#include <tbb/tbb.h>

template<typename UINT>
class TIBasis2DZ2
	: public Basis2D<UINT>
{
public:
	struct RepData
	{
		std::size_t idx;
		uint32_t numRep;
		int p;//Z2 parity
	};

private:
	const int parity_; //Z2_parity

	tbb::concurrent_vector<UINT> rpts_; //Representatives
	tbb::concurrent_unordered_map<UINT, RepData> repDatas_;

	void constructBasis(bool useU1)
	{
		tbb::concurrent_vector<std::tuple<UINT,uint32_t> > candids;

		if(useU1)
		{
			const unsigned int n = this->getN();
			const unsigned int nup = n/2;

			BasisJz<UINT> basis(n,nup);

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
				rpts_.emplace_back(rep);
				repDatas_[rep] = RepData{0, numRep, 0};
			}
			else if(flipedRep > rep)
			{
				rpts_.emplace_back(rep);
				repDatas_[rep] = RepData{0, numRep, 1};
			}
			else {
				;
			}
		});
		tbb::parallel_sort(rpts_.begin(), rpts_.end());
		tbb::parallel_for(static_cast<UINT>(0u), static_cast<UINT>(rpts_.size()), 
				[&](UINT idx){
			repDatas_[rpts_[idx]].idx = idx;
		});
	}


public:
	TIBasis2DZ2(uint32_t Lx, uint32_t Ly, uint32_t kx, uint32_t ky, bool useU1, int parity = 1)
		: Basis2D<UINT>(Lx, Ly, kx, ky), parity_(parity)
	{
		assert((parity == 1) || (parity == -1));
		constructBasis(useU1);
	}

	TIBasis2DZ2<UINT>(const TIBasis2DZ2<UINT>& ) = default;
	TIBasis2DZ2<UINT>(TIBasis2DZ2<UINT>&& ) = default;

	inline int getParity() const { return parity_; }

	std::size_t getDim() const override
	{
		return rpts_.size();
	}

	UINT getNthRep(int n) const override
	{
		return rpts_[n];
	}

	RepData getRepData(int n) const
	{
		return repDatas_[rpts_[n]];
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
		auto aRepData = repDatas_.at(rpts_[aidx]);
		double Na = aRepData.numRep/double(1+aRepData.p);

		UINT bRep;
		uint32_t bRotX, bRotY;
		std::tie(bRep, bRotX, bRotY) = this->getMinRots(bSigma);
		auto iter = repDatas_.find(bRep);
		if(iter == repDatas_.end())
		{
			c *= parity_;
			std::tie(bRep, bRotX, bRotY) = this->getMinRots(flip(bSigma));
			iter = repDatas_.find(bRep);
			if(iter == repDatas_.end())
				return std::make_pair(-1, 0.0);
		}

		auto bRepData = iter->second;
		double Nb = bRepData.numRep/double(1+bRepData.p);

		return std::make_pair(bRepData.idx, 
				c*sqrt(Nb/Na)*this->phase(bRotX,bRotY));
	}

	std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const override
	{
		using std::pow;

		std::map<UINT, double> res;
		UINT rep = getNthRep(n);
		auto repData = repDatas_.at(rep);
	
		double norm = 1.0/sqrt(repData.numRep*(1 + repData.p)*this->getN());

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
		if(repData.p != 0)
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
