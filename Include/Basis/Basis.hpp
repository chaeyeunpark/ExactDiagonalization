#ifndef ED_BASIS_HPP
#define ED_BASIS_HPP
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>
#include <initializer_list>


template<typename UINT>
class Basis
{
private:
	const unsigned int N_;
	const UINT ups_;
public:

	Basis(unsigned int N)
		: N_{N}, ups_{(UINT(1) << N_) -1}
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

	UINT rotl(UINT value, unsigned int count) const
	{
		count %= N_;
		return ((value << count) & ups_) | (value >> (N_ - count));
	}

	std::pair<UINT, unsigned int> findMinRots(UINT sigma) const
	{
		UINT rep = sigma;
		unsigned int rot = 0u;
		for(unsigned int r = 1; r < N_; r++)
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

	virtual std::size_t getDim() const = 0;
	virtual UINT getNthRep(int n) const = 0;
	virtual std::pair<int, double> hamiltonianCoeff(UINT bsigma, int adix) const = 0;

	virtual std::vector<std::pair<UINT, double>> basisVec(unsigned int n) const = 0;

	virtual ~Basis(){}
};

#endif//ED_BASIS_HPP
