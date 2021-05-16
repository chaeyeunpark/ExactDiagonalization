#pragma once

#include <map>
#include <vector>
#include <cassert>
#include <cmath>

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_sort.h>

#include "Basis.hpp"

struct RepData 
{
	std::size_t rptIdx;
	int rot;
	int parity;
};

template<typename UINT>
class TIBasisZ2
	: public Basis<UINT>
{
private:
	const int k_;
	const int p_;

	tbb::concurrent_vector<UINT> rpts_; //rpts_ is NOT sorted
	tbb::concurrent_unordered_map<UINT, RepData > parity_;

	int checkState(UINT s) const
	{
		UINT sr = s;
		const auto N = this->getN();
		for(int r = 1; r <= N; r++)
		{
			sr = this->rotl(s, r);
			if(sr < s)
			{
				return -1; //s is not a representative
			}
			else if(sr == s)
			{
				if((k_ % (N/r)) != 0)
					return -1; //this representative is not allowed for this k
				return r;
			}
		}
		return -1;
	}
	
	inline bool checkParity(int rot) const
	{
		if (p_ == 1)
			return ((rot*k_ % this->getN()) == 0);
		else //if p_ == -1
			return ((rot*k_ % this->getN()) == this->getN()/2);
	}

	void constructBasis()
	{
		tbb::concurrent_vector<std::pair<UINT, int> > ss;
		const unsigned int N = this->getN();
		ss.reserve((1<<(N-3))/N);
		{
			UINT s = 0;
			int r = checkState(s);
			if(r > 0)
			{
				ss.emplace_back(s,r);
			}
		}

		tbb::parallel_for(static_cast<UINT>(1), (UINT(1)<<UINT(N)), static_cast<UINT>(2), 
				[&](UINT s)
		{
			int r = checkState(s);
			if(r > 0)
			{
				ss.emplace_back(s,r);
			}
		});


		tbb::parallel_for(static_cast<std::size_t>(0), ss.size(), 
					[&](std::size_t idx)
		{
			UINT rep = ss[idx].first;
			auto s = this->findMinRots(flip(rep)); //s is the representative for flip(rep)
			if(s.first == rep && checkParity(s.second))
			{
				rpts_.emplace_back(rep);
				parity_[rep] = RepData{0, ss[idx].second, 0};
			}
			else if(s.first > rep)
			{
				rpts_.emplace_back(rep);
				parity_[rep] = RepData{0, ss[idx].second, 1};
			}
			else //s.first < rep
			{
				;
			}
		});

		//sort to make it consistent over different instances
		tbb::parallel_sort(rpts_);
		for(std::size_t idx = 0; idx < rpts_.size(); ++idx)
		{
			parity_[rpts_[idx]].rptIdx = idx;
		}

		//parity_ and rpts_ constructed

	}

public:
	TIBasisZ2(unsigned int N, unsigned int k, int p)
		: Basis<UINT>{N}, k_(k), p_(p)
	{
		assert(k == 0 || ((k == N/2) && (N%2 == 0)));
		assert(p_ == 1 || p_ == -1);
		constructBasis();
	}

	inline UINT flip(UINT value) const
	{
		return ((this->getUps())^value);
	}

	inline int getK() const { return k_; }
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
		double expk = (k_==0)?1.0:-1.0;

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
	std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const
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
