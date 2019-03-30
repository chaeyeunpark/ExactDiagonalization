#ifndef CY_TI_BASIS_Z2_HPP
#define CY_TI_BASIS_Z2_HPP

#include <map>
#include <vector>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>

#include "Basis.hpp"

struct RepData 
{
	std::size_t index;
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

	std::vector<UINT> rpts_;
	std::map<UINT, RepData> parity_;

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
					return -1; //this representative is not allowed for this k
				return r;
			}
		}
		return -1;
	}
	
	bool checkParity(int rot) const
	{
		if (p_ == 1)
			return ((rot*k_ % this->getN()) == 0);
		else
			return ((rot*k_ % this->getN()) == this->getN()/2);
	}

	void constructBasis()
	{
		std::vector<std::pair<UINT, int> > ss;
		{
			UINT s = 0;
			int r = checkState(s);
			if(r > 0)
			{
				ss.emplace_back(s,r);
			}
		}
		for(UINT s = 1; s <= this->getUps(); s+=2)
		{
			int r = checkState(s);
			if(r > 0)
			{
				ss.emplace_back(s,r);
			}
		}
		//now ss are sorted candidates
		for(auto iter = ss.begin(); iter != ss.end(); ++iter)
		{
			UINT rep = iter->first;
			auto s = this->findRepresentative(flip(rep));
			if(s.first == rep && checkParity(s.second))
			{
				rpts_.emplace_back(rep);
				parity_[rep] = RepData{rpts_.size()-1, iter->second, 0};
			}
			else if(s.first > rep)
			{
				rpts_.emplace_back(rep);
				parity_[rep] = RepData{rpts_.size()-1, iter->second, 1};
			}
			else //s.first < rep
			{
				;
			}
		}

	}

public:
	TIBasisZ2(int N, int k, int p)
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

	const std::map<UINT, RepData>& getData() const
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
		std::tie(bRep, bRot) = this->findRepresentative(bSigma);
		auto iter = parity_.find(bRep);
		if(iter == parity_.end())
		{
			c *= p_;
			std::tie(bRep, bRot) = this->findRepresentative(this->flip(bSigma));
			iter = parity_.find(bRep);

			if(iter == parity_.end())
				return std::make_pair(-1, 0.0);
		}
		auto pb = iter->second;
		double Nb = 1.0/double(1 + abs(pb.parity))/pb.rot;

		assert(pb.index < rpts_.size());

		return std::make_pair(pb.index,sqrt(Nb/Na)*pow(expk, bRot)*c);
	}
	Eigen::MatrixXd basisMatrix() const
	{
		double expk = (k_==0)?1.0:-1.0;
		Eigen::MatrixXd res(rpts_.size(), 1<<(this->getN()));
		res.setZero();
		for(int i = 0; i < rpts_.size(); i++)
		{
			auto rep = getNthRep(i);
			auto p = parity_.at(rep);
			for(int k = 0; k < p.rot; k++)
			{
				res(i, this->rotl(rep,k)) = pow(expk,k);
			}
			if(p.parity == 0)
			{
				res.row(i).normalize();
				continue;
			}
			rep = this->flip(rep);
			for(int k = 0; k < p.rot; k++)
			{
				res(i, this->rotl(rep,k)) = p_*pow(expk,k);
			}
			res.row(i).normalize();
		}
		return res;
	}
	/*
	template<typename T>
	T overlap(int idx, const Eigen::Matrix<T,  Eigen::Dynamic, 1>& rhs) const
	{
		double expk = (k_==0)?1.0:-1.0;
		auto rpt = rpts_[idx]; //representative
		auto p = parity_.at(rpts_[idx]);
		double Na = 1.0/double(1 + abs(p.parity))/p.rot;
		if(p.parity == 0)
		{
			T res = 0.0;
			for(int i = 0; i < p.rot; i++)
			{
				res += rhs(rotl(rpt, i))*Na*pow(expk,i);
			}
			return res;
		}
		else
		{
			T res = 0.0;
			for(int i = 0; i < p.rot; i++)
			{
				res += rhs(rotl(rpt, i))*Na*pow(expk,i);
				res -= rhs(rotl(rpt, i))*Na*pow(expk,i);
			}
			return res;
		}
	}*/
};
#endif//CY_TI_BASIS_HPP
