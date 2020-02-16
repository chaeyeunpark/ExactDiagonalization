#ifndef CY_TI_BASIS_HPP
#define CY_TI_BASIS_HPP

#include <boost/serialization/vector.hpp>
#include <cmath>
#include <iostream>

#include <Eigen/Dense>

#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include "Basis/Basis.hpp"

template<typename UINT>
class TIBasis
	: public Basis<UINT>
{
private:
	const int k_;

	tbb::concurrent_vector<std::pair<UINT, int> > rpts_; //Representatives


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
			UINT sr = this->rotl(sigma, r);
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
		auto rptComp = 
			[](const std::pair<UINT, int>& a, const std::pair<UINT, int>& b) -> bool
		{
			return a.first < b.first;
		};
		{
			UINT s = 0;
			int r = checkState(s);
			if(r > 0)
			{
				rpts_.emplace_back(s, r);
			}
		}
		tbb::parallel_for(static_cast<UINT>(1), this->getUps(), static_cast<UINT>(2),
				[&](UINT s)
		{
			int r = checkState(s);
			if(r > 0)
			{
				rpts_.emplace_back(s,r);
			}
		});

		tbb::parallel_sort(rpts_.begin(), rpts_.end(), rptComp);
	}

public:
	TIBasis(int N, int k)
		: Basis<UINT>(N), k_(k)
	{
		constructBasis();
	}

	TIBasis(const TIBasis& ) = default;
	TIBasis(TIBasis&& ) = default;

	TIBasis& operator=(const TIBasis& ) = default;
	TIBasis& operator=(TIBasis&& ) = default;

	inline int getK() const { return k_; }

	inline int rotRpt(int n) const
	{
		return rpts_[n].second;
	}
	
	int stateIdx(UINT rep) const
	{
		auto rptComp = 
			[](const std::pair<UINT, int>& a, UINT v) -> bool
		{
			return a.first < v;
		};
		auto iter = std::lower_bound(rpts_.begin(), rpts_.end(), rep, rptComp);
		return distance(rpts_.begin(), iter);
	}

	const tbb::concurrent_vector<std::pair<UINT, int> >& getRepresentatives() const
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


	std::pair<int, double> hamiltonianCoeff(UINT bSigma, int aidx) const override
	{
		using std::sqrt;
		using std::pow;

		double expk = (k_==0)?1.0:-1.0;

		UINT bRep;
		int bRot;
		std::tie(bRep, bRot) = this->getMinRots(bSigma);

		int bidx = stateIdx(bRep);

		if(bidx >= getDim())
		{
			return std::make_pair(-1, 0.0);
		}

		double Na = 1.0/rpts_[aidx].second;
		double Nb = 1.0/rpts_[bidx].second;

		return std::make_pair(bidx, sqrt(Nb/Na)*pow(expk, bRot));
	}

	Eigen::VectorXd basisVec(unsigned int n) const
	{
		const double expk = (k_==0)?1.0:-1.0;
		Eigen::VectorXd res(1<<(this->getN()));
		res.setZero();

		auto rep = rpts_[n].first;
		for(int k = 0; k < rpts_[n].second; k++)
		{
			res( this->rotl(rep,k)) = pow(expk,k);
		}
		res /= sqrt(rpts_[n].second);
		return res;
	}

};
#endif//CY_TI_BASIS_HPP
