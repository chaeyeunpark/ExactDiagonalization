#pragma once
#include <climits>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>
#include <initializer_list>

namespace edlib
{
template<typename UINT>
class AbstractBasis
{
private:
	const unsigned int N_;
	const UINT ups_;

public:
	AbstractBasis(unsigned int N)
		: N_{N}, ups_{(~UINT(0)) >> (sizeof(UINT)*CHAR_BIT - N)}
	{
	}

	inline int getN() const
	{
		return N_;
	}

	inline UINT getUps() const
	{
		return ups_;
	}

	UINT mask(std::initializer_list<uint32_t> pos) const
	{
		UINT s = 0;
		for(unsigned int p: pos)
			s ^= (UINT(1) << p);
		return s;
	}
	virtual std::size_t getDim() const = 0;
	virtual UINT getNthRep(int n) const = 0;
	virtual std::pair<int, double> hamiltonianCoeff(UINT bsigma, int adix) const = 0;

	virtual std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const = 0;

	virtual ~AbstractBasis(){}
};
} // namespace edlib
