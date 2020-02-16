#ifndef CY_TI_BASIS_HPP
#define CY_TI_BASIS_HPP

#include <boost/serialization/vector.hpp>
#include <cmath>
#include <iostream>

#include "Basis.hpp"

template<typename UINT>
class TIBasis
	: public Basis<UINT>
{
private:
	const int k_;

	std::vector<UINT> rpts_; //Representatives
	std::vector<int> rotRpts_; //R value for each representative

	int checkState(UINT s) 
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
					return -1; //this representative is not allowed for k
				return r;
			}
		}
		return -1;
	}

	std::pair<UINT, int> getMinRots(UINT sigma) const
	{
		UINT rep = sigma;
		int rot = 0;
		for(int r = 1; r < this->getN(); r++)
		{
			UINT sr = rotl(sigma, r);
			if(sr < rep)
			{
				rep = sr;
				rot = r;
			}
		}
		return std::make_pair(rep, rot);
	}

	void constructBasis()
	{
		{
			UINT s = 0;
			int r = checkState(s);
			if(r > 0)
			{
				rpts_.push_back(s);
				rotRpts_.push_back(r);
			}
		}
		for(UINT s = 1; s <= this->getUps(); s+=2)
		{
			int r = checkState(s);
			if(r > 0)
			{
				rpts_.push_back(s);
				rotRpts_.push_back(r);
			}
		}
	}

public:
	TIBasis(int N, int k)
		: Basis<UINT>(N), k_(k)
	{
		constructBasis();
	}

	TIBasis(const TIBasis& ) = default;
	/* Default move constructor is not supported in intel c compiler */
	//TIBasis(TIBasis&& ) = default;

	TIBasis& operator=(const TIBasis& ) = default;
	/* Default move assignment is not supported in intel c compiler */
	//TIBasis& operator=(TIBasis&& ) = default;

	inline int getK() const { return k_; }

	inline int rotRpt(int n) const
	{
		return rotRpts_[n];
	}
	
	std::size_t stateIdx(UINT rep) const
	{
		auto iter = lower_bound(rpts_.begin(), rpts_.end(), rep);
		return distance(rpts_.begin(), iter);
	}

	std::vector<UINT> getRepresentatives() const
	{
		return rpts_;
	}

	std::size_t getDim() const override
	{
		return rpts_.size();
	}

	UINT getNthRep(std::size_t n) const override
	{
		return rpts_[n];
	}

	std::pair<int, double> hamiltonianCoeff(UINT bSigma, std::size_t aidx) const override
	{
		using std::sqrt;
		using std::pow;

		double expk = (k_==0)?1.0:-1.0;

		UINT bRep;
		int bRot;
		std::tie(bRep, bRot) = this->getMinRots(bSigma);

		std::size_t bidx = stateIdx(bRep);

		if(bidx >= getDim())
		{
			return std::make_pair(-1, 0.0);
		}

		double Na = 1.0/rotRpts_[aidx];
		double Nb = 1.0/rotRpts_[bidx];

		return std::make_pair(bidx, sqrt(Nb/Na)*pow(expk, bRot));
	}
};
#endif//CY_TI_BASIS_HPP
