#pragma once
#include <map>
#include <vector>
#include <cassert>
#include <cmath>

#include <tbb/tbb.h>

#include "AbstractBasis1D.hpp"
#include "BasisJz.hpp"

namespace edlib
{
template<typename UINT>
class BasisZ2
	: public AbstractBasis1D<UINT>
{
private:
	const int N_;
	const int p_;

	tbb::concurrent_vector<std::pair<UINT, int>> rpts_;

	void constructBasis()
	{
		const unsigned int N = this->getN();

		tbb::parallel_for(UINT(0), UINT(1) << N, [&](std::size_t n)
		{
			auto s = flip(n);
			if(s == n)
			{
				rpts_.emplace_back(n, 0);
			}
			else if(s)
			{
				rpts_.emplace_back(n, 1);
			}
			else //s.first < rep
			{
				;
			}
		});

		//sort to make it consistent over different instances
		tbb::parallel_sort(rpts_, [](auto t1, auto t2){
			return t1.first < t2.first;
		});
		//parity_ and rpts_ constructed
	}

public:
	BasisZ2(unsigned int N, int p)
		: AbstractBasis<UINT>{N}, p_(p)
	{
		assert(p_ == 1 || p_ == -1);
		constructBasis();
	}

	BasisZ2(const BasisZ2& ) = default;
	BasisZ2(BasisZ2&& ) = default;

	inline UINT flip(UINT value) const
	{
		return ((this->getUps())^value);
	}

	inline int getP() const { return p_; }

	std::size_t getDim() const override
	{
		return rpts_.size();
	}

	UINT getNthRep(int n) const override
	{
		return rpts_[n].first;
	}

	std::pair<int,double> hamiltonianCoeff(UINT bSigma, int aidx) const override
	{
		UINT rep = rpts_[aidx].first;
		UINT pa = rpts_[aidx].second;
		double Na = 1.0/double(1 + abs(pa));

		double c = 1.0;

		auto iter = rpts_.find(bSigma);
		if(iter == rpts_.end())
		{
			c *= p_;
			std::tie(bRep, bRot) = this->findMinRots(this->flip(bSigma));
			iter = parity_.find(bRep);

			if(iter == parity_.end())
				return std::make_pair(-1, 0.0);
		}
		auto pb = iter->second;
		double Nb = 1.0/double(1 + abs(pb.parity))/pb.rot;

		return std::make_pair(pb.rptIdx,
				sqrt(Nb/Na)*pow(expk, bRot)*c);
	}
	

	/// return a vector of index/value pairs
	std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const override
	{
		const double expk = (k_==0)?1.0:-1.0;
		std::vector<std::pair<UINT,double>> res;

		auto rep = getNthRep(n);
		auto p = parity_.at(rep);
		double norm;
		if(p.parity == 0)
		{
			 norm = 1.0/sqrt(p.rot);
		}
		else
		{
			norm = 1.0/sqrt(2.0*p.rot);
		}

		for(int k = 0; k < p.rot; k++)
		{
			res.emplace_back(this->rotl(rep,k), pow(expk,k)*norm);
		}
		if(p.parity == 0)
		{
			return res;
		}
		rep = this->flip(rep);
		for(int k = 0; k < p.rot; k++)
		{
			res.emplace_back(this->rotl(rep,k), p_*pow(expk,k)*norm);
		}
		return res;
	}
};
} // namespace edlib
