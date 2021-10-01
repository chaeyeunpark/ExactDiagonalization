#pragma once
#include <cassert>
#include "AbstractBasis.hpp"

#include <tbb/tbb.h>

namespace edlib
{
template <typename UINT>
class BasisFullZ2
	: public AbstractBasis<UINT>
{
private:
	int parity_;
public:
	BasisFullZ2(unsigned int N, int parity)
		: AbstractBasis<UINT>(N), parity_{parity}
	{
		assert(parity == -1 || parity == 1);
	}

	std::size_t getDim() const override
	{
		return std::size_t(1) << (this->getN() - 1);
	}

	UINT getNthRep(int n) const override
	{
		return n;
	}

	std::pair<int, double> hamiltonianCoeff(UINT bsigma, int aidx) const override
	{
		if(bsigma < flip(bsigma))
			return std::make_pair(int(bsigma), 1.0);
		else
			return std::make_pair(int(flip(bsigma)), parity_);
	}

	std::vector<std::pair<UINT,double> > basisVec(unsigned int n) const override
	{
		using std::sqrt;
		return std::vector<std::pair<UINT, double>>{
			std::make_pair(UINT(n), 1.0/sqrt(2.0)),
			std::make_pair(UINT(flip(n)), parity_*1.0/sqrt(2.0)),
		};
	}

	inline UINT flip(UINT value) const
	{
		return ((this->getUps())^value);
	}
};
} //namespace edlib
