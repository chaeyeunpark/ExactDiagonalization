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
class Basis1DZ2
	: public AbstractBasis1D<UINT>
{
public:
	struct RepData 
	{
		std::size_t rptIdx;
		int rot;
		int parity;
	};

private:
	const int p_;

	tbb::concurrent_vector<UINT> rpts_; 
	tbb::concurrent_unordered_map<UINT, RepData> parity_;

	int phase(int rot) const
	{
		const auto k = this->getK();
		int N = this->getN();
		if ((k*rot % N) == 0)
			return 1;
		else
			return -1;
	}

	void constructBasis(bool useU1)
	{
		tbb::concurrent_vector<std::pair<UINT, int> > candids;
		const unsigned int N = this->getN();
		candids.reserve((1<<(N-3))/N);

		if(useU1)
		{
			const unsigned int nup = N/2;

			BasisJz<UINT> basis(N, nup);

			tbb::parallel_for_each(basis.begin(), basis.end(), [&](UINT s)
			{
				int r = this->checkState(s);
				if(r > 0)
				{
					candids.emplace_back(s, r);
				}
			});
		}
		else
		{
			{//insert 0
				UINT s = 0;
				int r = this->checkState(s);
				if(r > 0)
				{
					candids.emplace_back(s, r);
				}
			}

			//insert other candids. Exclude even as their LSB is 0 (shifted one must be smaller)
			tbb::parallel_for(static_cast<UINT>(1), (UINT(1)<<UINT(N)), static_cast<UINT>(2), 
					[&](UINT s)
			{
				int r = this->checkState(s);
				if(r > 0)
				{
					candids.emplace_back(s, r);
				}
			});
		}

		tbb::parallel_for(static_cast<std::size_t>(0), candids.size(), 
					[&](std::size_t idx)
		{
			UINT rep = candids[idx].first;
			auto s = this->findMinRots(flip(rep));
			if(s.first == rep && phase(s.second)*p_ == 1)
			{
				rpts_.emplace_back(rep);
				parity_.emplace(rep, RepData{0, candids[idx].second, 0});
			}
			else if(s.first > rep)
			{
				rpts_.emplace_back(rep);
				parity_.emplace(rep, RepData{0, candids[idx].second, 1});
			}
			else //s.first < rep
			{
				;
			}
		});

		//sort to make it consistent over different instances
		tbb::parallel_sort(rpts_);
		tbb::parallel_for(static_cast<UINT>(0u), static_cast<UINT>(rpts_.size()), 
				[&](UINT idx){
			parity_[rpts_[idx]].rptIdx = idx;
		});


		//parity_ and rpts_ constructed

	}

public:
	Basis1DZ2(unsigned int N, unsigned int k, int p, bool useU1)
		: AbstractBasis1D<UINT>{N, k}, p_(p)
	{
		assert(k == 0 || ((k == N/2) && (N%2 == 0)));
		assert(p_ == 1 || p_ == -1);
		constructBasis(useU1);
	}

	Basis1DZ2(const Basis1DZ2& ) = default;
	Basis1DZ2(Basis1DZ2&& ) = default;

	inline UINT flip(UINT value) const
	{
		return ((this->getUps())^value);
	}

	inline int getP() const { return p_; }

	RepData getData(UINT s) const
	{
		return parity_.at(s);
	}

	const tbb::concurrent_unordered_map<UINT, RepData >& getParityMap() const
	{
		return parity_;
	}

	std::size_t getDim() const override
	{
		return rpts_.size();
	}

	UINT getNthRep(int n) const override
	{
		return rpts_[n];
	}

	std::pair<int,double> hamiltonianCoeff(UINT bSigma, int aidx) const override
	{
		const auto k = this->getK();
		double expk = (k == 0)?1.0:-1.0;

		auto pa = parity_.at(rpts_[aidx]);
		double Na = 1.0/double(1 + abs(pa.parity))/pa.rot;

		double c = 1.0;

		UINT bRep;
		int bRot;
		std::tie(bRep, bRot) = this->findMinRots(bSigma);
		auto iter = parity_.find(bRep);
		if(iter == parity_.end())
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
		const auto k = this->getK();
		const double expk = (k == 0)?1.0:-1.0;
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