#pragma once
#include <boost/serialization/vector.hpp>
#include <cmath>
#include <iostream>

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_for_each.h>

#include "BasisJz.hpp"
#include "AbstractBasis1D.hpp"

namespace edlib
{
template<typename UINT>
class Basis1D
	: public AbstractBasis1D<UINT>
{
private:
	tbb::concurrent_vector<std::pair<UINT, unsigned int>> rpts_; //Representatives

	int checkState(UINT s) 
	{
		UINT sr = s;
		const auto N = this->getN();
		const auto k = this->getK();

		for(int r = 1; r <= N; r++)
		{
			sr = this->rotl(s, r);
			if(sr < s)
			{
				return -1; //s is not a representative
			}
			else if(sr == s)
			{
				/* s is smller than rotl(s,1), rot(s, 2), ..., rot(s, r-1)
				 * when we fall in this else if clause. As rot(s,r) == s,
				 * s is the smallest among rotl(s, 1), ..., rotl(s, N-1).
				 */
				if((k % (N/r)) != 0)
					return -1; //this representative is not allowed for k
				return r;
			}
		}
		return -1;
	}

	void constructBasisFull()
	{
		// insert 0
		{
			UINT s = 0;
			int r = checkState(s);
			if(r > 0)
			{
				rpts_.emplace_back(s,r);
			}
		}

		// iterate over all odd numbers
		const unsigned int N = this->getN();
		tbb::parallel_for(UINT(1), (UINT(1)<<UINT(N)), UINT(2), [&](UINT s)
		{
			int r = checkState(s);
			if(r > 0)
			{
				rpts_.emplace_back(s,r);
			}
		});
		tbb::parallel_sort(rpts_.begin(), rpts_.end());
	}

	void constructBasisJz()
	{
		const unsigned int n = this->getN();
		const unsigned int nup = n/2;

		BasisJz<UINT> basis(n,nup);

		tbb::parallel_for_each(basis.begin(), basis.end(), [&](UINT s)
		{
			int r = checkState(s);
			if(r > 0)
			{
				rpts_.emplace_back(s,r);
			}
		});
		tbb::parallel_sort(rpts_.begin(), rpts_.end());
	}

public:
	Basis1D(unsigned int N, unsigned int k, bool useU1)
		: AbstractBasis1D<UINT>(N, k)
	{
		if(useU1)
		{
			constructBasisJz();
		}
		else
		{
			constructBasisFull();
		}
	}

	Basis1D(const Basis1D& ) = default;
	Basis1D(Basis1D&& ) = default;

	//Basis1D& operator=(const Basis1D& ) = default;
	//Basis1D& operator=(Basis1D&& ) = default;

	unsigned int stateIdx(UINT rep) const
	{
		auto comp = [](const std::pair<UINT, unsigned int>& v1, UINT v2)
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

	tbb::concurrent_vector<std::pair<UINT,unsigned int> > getRepresentatives() const
	{
		return rpts_;
	}

	std::size_t getDim() const override
	{
		return rpts_.size();
	}

	UINT getNthRep(int n) const override
	{
		return rpts_[n].first;
	}

	inline unsigned int rotRpt(int n) const
	{
		return rpts_[n].second;
	}

	std::pair<int, double> hamiltonianCoeff(UINT bSigma, int aidx) const override
	{
		using std::sqrt;
		using std::pow;

		const auto k = this->getK();
		double expk = (k == 0)?1.0:-1.0;

		UINT bRep;
		int bRot;
		std::tie(bRep, bRot) = this->findMinRots(bSigma);

		auto bidx = stateIdx(bRep);

		if(bidx >= getDim())
		{
			return std::make_pair(-1, 0.0);
		}

		double Na = 1.0/rpts_[aidx].second;
		double Nb = 1.0/rpts_[bidx].second;

		return std::make_pair(bidx, sqrt(Nb/Na)*pow(expk, bRot));
	}

	std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const
	{
		const auto k = this->getK();
		const double expk = (k == 0)?1.0:-1.0;
		std::vector<std::pair<UINT,double>> res;

		auto rep = rpts_[n].first;
		double norm = 1.0/sqrt(rpts_[n].second);
		for(unsigned int r = 0; r < rpts_[n].second; r++)
		{
			res.emplace_back( this->rotl(rep, r), pow(expk, r)*norm);
		}
		return res;
	}
};
} // namespace edlib
