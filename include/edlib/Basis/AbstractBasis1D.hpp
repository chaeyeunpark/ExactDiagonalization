#pragma once
#include "AbstractBasis.hpp"

namespace edlib
{
template<typename UINT>
class AbstractBasis1D
	: public AbstractBasis<UINT>
{
protected:
	const unsigned int k_;

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

public:
	AbstractBasis1D(unsigned int N, unsigned int k)
		: AbstractBasis<UINT>(N), k_{k}
	{
	}

	UINT rotl(UINT value, unsigned int count) const
	{
		const auto N = this->getN();
		count %= N;
		return ((value << count) & this->getUps()) | (value >> (N - count));
	}

	std::pair<UINT, unsigned int> findMinRots(UINT sigma) const
	{
		const auto N = this->getN();
		UINT rep = sigma;
		unsigned int rot = 0u;
		for(unsigned int r = 1; r < N; r++)
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

	inline unsigned int getK() const { return k_; }
};
} // namespace edlib
