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
	const uint32_t N_;
	const UINT ups_;

public:
	AbstractBasis(uint32_t N)
		: N_{N}, ups_{(~UINT(0)) >> (sizeof(UINT)*CHAR_BIT - N)}
	{
	}

	inline uint32_t getN() const
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
		for(uint32_t p: pos)
			s ^= (UINT(1) << p);
		return s;
	}

	virtual std::size_t getDim() const = 0;
	virtual UINT getNthRep(uint32_t n) const = 0;
	virtual std::pair<int, double> hamiltonianCoeff(UINT bsigma, int adix) const = 0;

	virtual std::vector<std::pair<UINT, double>> basisVec(uint32_t n) const = 0;

	virtual ~AbstractBasis() = default;
};
} // namespace edlib
