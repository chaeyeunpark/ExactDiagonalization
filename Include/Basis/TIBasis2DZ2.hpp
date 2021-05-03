#pragma once

#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cassert>

#include "Basis.hpp"
#include "BasisJz.hpp"

#include <tbb/tbb.h>

template<typename UINT>
class TIBasis2DZ2
	: public Basis<UINT>
{
public:
	struct RepData
	{
		std::size_t idx;
		uint32_t numRep;
		int p;//Z2 parity
	};

private:
	const uint32_t Lx_;
	const uint32_t Ly_;
	const uint32_t kx_;
	const uint32_t ky_;
	const int parity_; //Z2_parity
	const UINT px_;

	tbb::concurrent_vector<UINT> rpts_; //Representatives
	tbb::concurrent_unordered_map<UINT, RepData> repDatas_;

	/**
	 * In two dimension, there are several different rotations that gives 
	 * the same representation
	 * */
	std::tuple<UINT, uint32_t, uint32_t> getMinRots(UINT sigma) const
	{
		UINT rep = sigma;
		uint32_t rotX = Lx_;
		uint32_t rotY = Ly_;
		for(uint32_t rx = 1; rx <= Lx_; rx++)
		{
			UINT srx = rotateX(sigma, rx);
			for(uint32_t ry = 1; ry <= Ly_; ry++)
			{
				UINT sr = rotateY(srx, ry);
				if(sr < rep)
				{
					rep = sr;
					rotX = rx;
					rotY = ry;
				}
			}
		}
		return std::make_tuple(rep, rotX, rotY);
	}

	uint32_t checkState(UINT s)  const
	{
		uint32_t cnt = 0;
		UINT sr = s;

		for(uint32_t rx = 1; rx <= Lx_; rx++)
		{
			for(uint32_t ry = 1; ry <= Ly_; ry++)
			{
				sr = rotateX(s, rx);
				sr = rotateY(sr, ry);
				if(sr < s)
				{
					return 0; //s is not a representative
				}
				else if(sr == s)
				{
					if(phase(rx,ry) == -1) //not allowed
					{
						return 0;
					}
					else
					{
						++cnt;
					}
				}
			}
		}
		return cnt;
	}

	int phase(int32_t rotX, int32_t rotY) const
	{
		//return exp(2\Pi I(kx_*rotX/Lx_ + ky_*rotY/Ly_))
		int sgn = 1;
		if((kx_*rotX) % Lx_ != 0)
		{
			sgn *= -1;
		}
		if((ky_*rotY) % Ly_ != 0)
		{
			sgn *= -1;
		}
		return sgn;
	}

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
				uint32_t numRep = checkState(s);
				if(numRep > 0)
				{
					candids.emplace_back(s, numRep);
				}
			});
		}
		else
		{
			//if the MSB is 1, the fliped one is smaller.
			tbb::parallel_for(UINT(0u), (UINT(1u) << UINT(Lx_*Ly_ - 1)),
			[&](UINT s) {
				uint32_t numRep = checkState(s);
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
			std::tie(flipedRep, rotX, rotY) = getMinRots(fliped);

			if(flipedRep == rep && phase(-rotX, -rotY)*parity_ == 1)
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
		: Basis<UINT>(Lx*Ly), Lx_(Lx), Ly_(Ly),
		kx_(kx), ky_(ky), parity_(parity),
		px_{(~UINT(0)) >> (sizeof(UINT)*8 - Lx_)}
	{
		assert((parity == 1) || (parity == -1));
		constructBasis(useU1);
	}

	TIBasis2DZ2<UINT>(const TIBasis2DZ2<UINT>& ) = default;
	TIBasis2DZ2<UINT>(TIBasis2DZ2<UINT>&& ) = default;

	UINT rotateY(UINT sigma, int32_t r) const
	{
		return this->rotl(sigma, Lx_*r);
	}

	UINT rotateX(UINT sigma, int32_t r) const
	{
		assert((r >= 0) && (r <= Lx_));
		const auto Lx = Lx_;
		const auto Ly = Ly_;
		const auto px = px_;
		for(uint32_t ny = 0; ny < Ly; ny++)
		{
			UINT t = (sigma >> (ny*Lx)) & px;
			t = ((t << r) | (t >> (Lx-r))) & px;
			sigma &= ~(px << ny*Lx);
			sigma |= (t << ny*Lx);
		}
		return sigma;
	}





	inline uint32_t getLx() const { return Lx_; }
	inline uint32_t getLy() const { return Ly_; }
	inline uint32_t getKx() const { return kx_; }
	inline uint32_t getKy() const { return ky_; }
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

	inline uint32_t toIdx(uint32_t nx, uint32_t ny) const 
	{
		return ny*Lx_ + nx;
	}

	std::pair<int, double> hamiltonianCoeff(UINT bSigma, int aidx) const override
	{
		using std::abs;
		using std::sqrt;
		using std::pow;

		double c = 1.0;
		auto aRepData = repDatas_.at(rpts_[aidx]);
		double Na = 1.0/double(1+aRepData.p)*aRepData.numRep;

		UINT bRep;
		uint32_t bRotX, bRotY;
		std::tie(bRep, bRotX, bRotY) = getMinRots(bSigma);
		auto iter = repDatas_.find(bRep);
		if(iter == repDatas_.end())
		{
			c *= parity_;
			std::tie(bRep, bRotX, bRotY) = getMinRots(flip(bSigma));
			iter = repDatas_.find(bRep);
			if(iter == repDatas_.end())
				return std::make_pair(-1, 0.0);
		}

		auto bRepData = iter->second;
		double Nb = 1.0/double(1+bRepData.p)*bRepData.numRep;

		return std::make_pair(bRepData.idx, 
				c*sqrt(Nb/Na)*phase(bRotX,bRotY));
	}

	std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const override
	{
		using std::pow;

		std::map<UINT, double> res;
		UINT rep = getNthRep(n);
		auto repData = repDatas_.at(rep);
	
		double norm = 1.0/sqrt(repData.numRep*(1 + repData.p)*this->getN());

		for(uint32_t nx = 0; nx < Lx_; nx++)
		{
			auto srx = rotateX(rep, nx);
			for(uint32_t ny = 0; ny < Ly_; ny++)
			{
				auto sr = rotateY(srx, ny);
				res[sr] += phase(nx, ny)*norm;
			}
		}
		if(repData.p != 0)
		{
			UINT fliped = flip(rep);
			for(uint32_t nx = 0; nx < Lx_; nx++)
			{
				auto srx = rotateX(fliped, nx);
				for(uint32_t ny = 0; ny < Ly_; ny++)
				{
					auto sr = rotateY(srx, ny);
					res[sr] += phase(nx, ny)*parity_*norm;
				}
			}
		}

		return std::vector<std::pair<UINT, double>>(res.begin(), res.end());
	}
};
