#pragma once
#include "AbstractBasis.hpp"
#include "BasisJz.hpp"

#include <tbb/tbb.h>

namespace edlib
{
template <typename UINT>
class FullBasis
	: AbstractBasis<UINT>
{
private:

public:
	FullBasis(unsigned int N)
		: AbstractBasis<UINT>(N)
	{
	}

	std::size_t getDim() const override
	{
		return (1u << this->getN());
	}

	UINT getNthRep(int n) const override
	{
		return UINT(n);
	}

	std::pair<int, double> hamiltonianCoeff(UINT bsigma, int aidx) const override
	{
		return std::make_pair(int(bsigma), 1.0);
	}

	std::vector<std::pair<UINT,double> > basisVec(unsigned int n) const override
	{
		return std::vector<std::pair<UINT,double>>{std::make_pair(UINT(n), 1.0)};
	}
};
} //namespace edlib
