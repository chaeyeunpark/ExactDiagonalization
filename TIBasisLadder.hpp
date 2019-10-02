#ifndef CY_TIBASISLADDER_HPP
#define CY_TIBASISLADDER_HPP

#include <boost/serialization/vector.hpp>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>

#include "Basis.hpp"

#ifndef NDEBUG
#include <Eigen/Dense>
#endif

namespace ED
{

template<typename UINT>
class TIBasisLadder
	: public Basis<UINT>
{
public:
	struct RepData
	{
		int rot;
		std::vector<uint32_t> qs; 
		/* each q consists of 3 bits. 
		 * q & 1 => parity_
		 * (q >> 1) & 1 => px;
		 * (q >> 2) & 1 => py;
		 * Acting the operators indicated by the set bits give the same translational minimum
		 * e.g. if s = 0101, px and py&parity_ preserve the representation => 
		 * qs = {0b000, 0b010, 0b101, 0b111}
		 */
	};

	struct ToRep
	{
		std::size_t idx; //index of real representative
		uint32_t p;
		/* p presents which operation should be taken to make it the rotation of real rep
		 */
	};

private:
	const int k_;
	const int parity_;
	const int px_;
	const int py_;

	std::map<UINT, ToRep> data_;

	std::vector<std::pair<UINT, RepData> > rpts_; //Representatives


	std::pair<UINT, int> getMinRots(UINT sigma) const
	{
		UINT rep = sigma;
		int rot = 0;
		for(int r = 1; r < this->getN()/2; r++)
		{
			UINT sr = this->rotl(sigma, 2*r);
			if(sr < rep)
			{
				rep = sr;
				rot = r;
			}
		}
		return std::make_pair(rep, rot);
	}


	int checkState(UINT s) 
	{
		UINT sr = s;
		const auto N = this->getN();
		const auto L = N/2;
		for(int r = 1; r <= L; r++)
		{
			sr = this->rotl(s, 2*r);
			if(sr < s)
			{
				return -1; //s is not a representative
			}
			else if(sr == s)
			{
				if((k_ % (L/r)) != 0)
					return -1; //this representative is not allowed for k
				return r;
			}
		}
		return -1;
	}

	void constructBasis()
	{
		std::vector<std::pair<UINT,int> > candids;
		{
			UINT s = 0u;
			int r = checkState(s);
			if(r > 0)
			{
				candids.emplace_back(s,r);
			}
		}

#pragma omp parallel
		{
#pragma omp for schedule(dynamic,8)
			for(UINT s = 1u; s <= this->getUps(); s+=4)
			{
				int r = checkState(s);
				if(r > 0)
				{
#pragma omp critical
					{
						candids.emplace_back(s,r);
					}
				}
			}
#pragma omp for schedule(dynamic,8)
			for(UINT s = 2u; s <= this->getUps(); s+=4)
			{
				int r = checkState(s);
				if(r > 0)
				{
#pragma omp critical
					{
						candids.emplace_back(s,r);
					}
				}
			}
#pragma omp for schedule(dynamic,8)
			for(UINT s = 3u; s <= this->getUps(); s+=4)
			{
				int r = checkState(s);
				if(r > 0)
				{
#pragma omp critical
					{
						candids.emplace_back(s,r);
					}
				}
			}
		}
		std::sort(candids.begin(), candids.end());

		//candids are sorted
		for(auto s: candids)
		{
			if(data_.find(s.first) != data_.end())
				continue;

			auto ss = reflections(s.first);

			std::array<std::pair<UINT, int>, 7> tiReps;
			
			std::transform(ss.begin()+1, ss.end(), tiReps.begin(), [this](UINT e)
			{
				return this->getMinRots(e);
			});

			std::array<bool, 7> leq;
			std::transform(tiReps.begin(), tiReps.end(), leq.begin(), [&s](const std::pair<UINT, int>& e) -> bool
			{
				return s.first <= e.first;
			});

			int L = this->getN()/2;
			if(std::all_of(leq.begin(), leq.end(), [](bool x){return x;}))
			{
				//s is true representative (smallest among all possible reflections)

				//check represntative is allowed for given quantum numbers
				bool isAllowed = true;
				std::vector<uint32_t> qs;
				qs.emplace_back(0);
				for(int i = 1; i < 8; i++)
				{
					if(tiReps[i-1].first == s.first)
					{
						if(((k_*tiReps[i-1].second % L) == L/2) && (getP(i) == 1))
						{
							isAllowed = false;
							break;
						}
						else if(((k_*tiReps[i-1].second % L) == 0) && (getP(i) == -1))
						{
							isAllowed = false;
							break;
						}
						else
						{
							qs.emplace_back(i);
						}
					}
				}
				if(!isAllowed)
					continue;

				auto idx = rpts_.size();
				rpts_.emplace_back(std::make_pair(s.first, RepData{s.second, std::move(qs)}));
				data_.emplace(std::make_pair(s.first, ToRep{idx, 0}));

				int bit_ordered[7] = {1, 2, 4, 3, 5, 6, 7};
				for(std::size_t i = 0; i < tiReps.size(); i++)
				{
					//if there are several tiReps[i] with the same value,
					//p with lowest bitcount will be inserted
					data_.emplace(std::make_pair(tiReps[bit_ordered[i]-1].first, ToRep{idx,bit_ordered[i]}));
				}
			}
		}
	}



public:
	TIBasisLadder(int N, int k, int parity, int px, int py)
		: Basis<UINT>(N), k_(k), parity_(parity), px_(px), py_(py)
	{
		assert(N%2 == 0);
		assert(parity == 1 || parity == -1);
		assert(px == 1 || px == -1);
		assert(py == 1 || py == -1);
		assert(k < N/2 && k >= 0);

		constructBasis();
	}

	TIBasisLadder(const TIBasisLadder& ) = default;
	TIBasisLadder(TIBasisLadder&& ) = default;

	TIBasisLadder& operator=(const TIBasisLadder& ) = default;
	TIBasisLadder& operator=(TIBasisLadder&& ) = default;

	inline int getK() const { return k_; }
	inline int getParity() const { return parity_; }
	inline int getPx() const { return px_; }
	inline int getPy() const { return py_; }

	int getP(uint32_t p) const
	{
		int P = 1;
		if(p & 1)
			P *= parity_;
		if( (p >> 1) & 1)
			P *= px_;
		if( (P >> 2) & 1)
			P *= py_;
		return P;
	}



	std::size_t getDim() const override
	{
		return rpts_.size();
	}

	UINT getNthRep(int n) const override
	{
		return rpts_[n].first;
	}

	RepData getRepData(int n) const
	{
		return rpts_[n].second;
	}

	inline UINT flip(UINT value) const
	{
		return ((this->getUps())^value);
	}

	inline UINT reflectX(UINT s) const
	{
		const auto N = this->getN();
		const auto L = N/2;
		UINT res = 0u;
		for(int i = 0; i < L; i++)
		{
			res |= (((s >> 2*i) & 0x3) << 2*(L-i-1));
		}
		return res;
	}



	inline UINT reflectY(UINT s) const
	{
		const auto N = this->getN();
		const auto L = N/2;
		constexpr UINT bt[] = {0,2,1,3};
		for(int i = 0; i < L; i++)
		{
			const UINT t = bt[(s >> 2*i) & 3];
			s = (s & ~(3UL << 2*i)) | (t << 2*i);
		}
		return s;
	}

	std::array<UINT,8> reflections(UINT s) const
	{
		auto s1 = flip(s);

		auto s2 = reflectX(s);
		auto s3 = reflectX(s1);

		auto s4 = reflectY(s);
		auto s5 = reflectY(s1);

		auto s6 = reflectY(s2);
		auto s7 = reflectY(s3);
		return {s,s1,s2,s3,s4,s5,s6,s7};
	}

	UINT reflect(UINT s, uint32_t p) const
	{
		switch(p & 1)
		{
			case 0:
				break;
			case 1:
				s = flip(s);
				break;
		}
		switch( (p >> 1) & 1)
		{
			case 0:
				break;
			case 1:
				s = reflectX(s);
				break;
		}
		switch( (p >> 2) & 1)
		{
			case 0:
				break;
			case 1:
				s = reflectY(s);
				break;
		}
		return s;
	}

	std::pair<int, double> hamiltonianCoeff(UINT bSigma, int aidx) const override
	{
		using std::sqrt;
		using std::pow;

		double expk = (k_==0)?1.0:-1.0;

		UINT bRep;
		int bRot;
		std::tie(bRep, bRot) = getMinRots(bSigma);

		auto iter = data_.find(bRep);
		if(iter == data_.end())
		{
			return std::make_pair(-1, 0.0);
		}

		int bidx = (iter->second).idx;
		UINT s = reflect(bSigma, (iter->second).p);
		std::tie(bRep, bRot) = getMinRots(s);

		assert(bRep == rpts_[bidx].first);

		double coeff = getP((iter->second).p)*pow(expk,bRot);

		double Na = sqrt(rpts_[aidx].second.qs.size())/sqrt(rpts_[aidx].second.rot);
		double Nb = sqrt(rpts_[bidx].second.qs.size())/sqrt(rpts_[bidx].second.rot);

		return std::make_pair(bidx, coeff*Nb/Na);
	}

#ifndef NDEBUG
	std::pair<int, double> coeff(UINT bSigma) const
	{
		using std::sqrt;
		using std::pow;

		const int N = this->getN();
		const int L = this->getN()/2;

		double expk = (k_==0)?1.0:-1.0;

		UINT bRep;
		int bRot;
		std::tie(bRep, bRot) = getMinRots(bSigma);

		auto iter = data_.find(bRep);
		if(iter == data_.end())
		{
			return std::make_pair(-1, 0.0);
		}

		int bidx = (iter->second).idx;
		UINT s = reflect(bSigma, (iter->second).p);
		std::tie(bRep, bRot) = getMinRots(s);

		assert(bRep == rpts_[bidx].first);

		double coeff = getP((iter->second).p)*pow(expk,bRot);

		double Nb = sqrt(rpts_[bidx].second.rot*8.0/rpts_[bidx].second.qs.size());

		return std::make_pair(bidx, coeff/Nb);
	}
#endif
};
}//namespace ED
#endif//CY_TIBASISLADDER_HPP
