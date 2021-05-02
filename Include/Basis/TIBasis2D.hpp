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
class TIBasis2D
	: public Basis<UINT>
{
private:
	const uint32_t Lx_;
	const uint32_t Ly_;
	const uint32_t kx_;
	const uint32_t ky_;
	const UINT px_;


	tbb::concurrent_vector< std::pair<UINT, uint32_t> > rpts_; //Representatives and number of repretations

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
					rpts_.emplace_back(s, numRep);
				}
			});
		}
		else
		{
			//if the MSB is 1, the fliped one is smaller.
			tbb::parallel_for(UINT(0u), (UINT(1u) << UINT(Lx_*Ly_)),
			[&](UINT s) {
				uint32_t numRep = checkState(s);
				if(numRep > 0)
				{
					rpts_.emplace_back(s, numRep);
				}
			});
		}

		tbb::parallel_sort(rpts_.begin(), rpts_.end());
	}


public:
	TIBasis2D(uint32_t Lx, uint32_t Ly, uint32_t kx, uint32_t ky, bool useU1)
		: Basis<UINT>(Lx*Ly), Lx_(Lx), Ly_(Ly),
		kx_(kx), ky_(ky),
		px_{(~UINT(0)) >> (sizeof(UINT)*8 - Lx_)}
	{
		constructBasis(useU1);
	}

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


	TIBasis2D<UINT>(const TIBasis2D<UINT>& ) = default;

	TIBasis2D<UINT>& operator=(const TIBasis2D<UINT>& ) = default;

	inline uint32_t getLx() const { return Lx_; }
	inline uint32_t getLy() const { return Ly_; }
	inline uint32_t getKx() const { return kx_; }
	inline uint32_t getKy() const { return ky_; }

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

	inline uint32_t toIdx(uint32_t nx, uint32_t ny) const 
	{
		return ny*Lx_ + nx;
	}

	std::pair<int, double> hamiltonianCoeff(UINT bSigma, int aidx) const override
	{
		using std::abs;
		using std::sqrt;
		using std::pow;

		uint32_t aNumRep = rpts_[aidx].second;
		double Na = aNumRep;

		UINT bRep;
		uint32_t bRotX, bRotY;
		std::tie(bRep, bRotX, bRotY) = getMinRots(bSigma);
		auto iter = std::lower_bound(rpts_.begin(), rpts_.end(), std::make_pair(bRep, 0u));
		if(iter == rpts_.end() || iter->first != bRep)
		{
			return std::make_pair(-1, 0.0);
		}
		
		int idx = std::distance(rpts_.begin(), iter);
		uint32_t bNumRep = iter->second;
		double Nb = bNumRep;

		return std::make_pair(idx, sqrt(Nb/Na)*phase(bRotX,bRotY));
	}

	std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const override
	{
		using std::pow;

		std::map<UINT, double> res;
		UINT rep = rpts_[n].first;
		uint32_t numRep = rpts_[n].second;
	
		double norm = 1.0/sqrt(numRep*this->getN());

		for(uint32_t nx = 0; nx < Lx_; nx++)
		{
			auto srx = rotateX(rep, nx);
			for(uint32_t ny = 0; ny < Ly_; ny++)
			{
				auto sr = rotateY(srx, ny);
				res[sr] += phase(nx, ny)*norm;
			}
		}

		return std::vector<std::pair<UINT, double>>(res.begin(), res.end());
	}

};
