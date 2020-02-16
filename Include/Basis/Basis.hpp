#ifndef ED_BASIS_HPP
#define ED_BASIS_HPP
#include <cstddef>
#include <cstdint>
#include <tuple>

template<typename UINT>
class Basis
{
private:
	const int N_;
	const UINT ups_;
public:

	Basis(int N) 
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

	UINT rotl(UINT value, unsigned int count) const
	{
		count %= N_;
		return ((value << count) & ups_) | (value >> (N_ - count));
	}

	std::pair<UINT, int> findRepresentative(UINT sigma) const
	{
		UINT rep = sigma;
		int rot = 0;
		for(int r = 1; r < N_; r++)
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
	virtual UINT getNthRep(std::size_t n) const = 0;
	/**
	 * \sum_i h_i|rep[aidx]> = coeff*|rep of bsigma>
	 */
	virtual std::pair<int, double> hamiltonianCoeff(UINT bsigma, std::size_t adix) const = 0;

	virtual ~Basis(){}
};

#endif//ED_BASIS_HPP
