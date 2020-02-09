#ifndef CY_TIXXXJ1J2_HPP
#define CY_TIXXXJ1J2_HPP
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <map>
#include <boost/dynamic_bitset.hpp>
//#include "BitOperations.h"
#include "Basis/Basis.hpp"

template<typename UINT>
class TIXXXJ1J2
{
private:
	Basis<UINT>& basis_;
	double J1_;
	double J2_;

public:
	TIXXXJ1J2(Basis<UINT>& basis, double J1, double J2)
		: basis_(basis), J1_(J1), J2_(J2)
	{
		
	}

	std::map<int,double> getCol(UINT n) const
	{
		int N = basis_.getN();

		UINT a = basis_.getNthRep(n);
		const boost::dynamic_bitset<> bs(N, a);

		std::map<int, double> m;
		for(int i = 0; i < N; i++)
		{
			//Next-nearest
			{
				int j = (i+1)%N;
				int sgn = (1-2*bs[i])*(1-2*bs[j]);

				m[n] += J1_*sgn;
				
				UINT s = a;
				UINT t = (1<<i) | (1<<j);
				s ^= t;

				int bidx;
				double coeff;

				std::tie(bidx, coeff) = basis_.hamiltonianCoeff(s, n);
				
				if(bidx >= 0)
					m[bidx] += J1_*(1-sgn)*coeff;
			}
			//Next-next-nearest
			{
				int j = (i+2)%N;
				int sgn = (1-2*bs[i])*(1-2*bs[j]);

				m[n] += J2_*sgn;
				
				UINT s = a;
				UINT t = (1<<i) | (1<<j);
				s ^= t;
		
				int bidx;
				double coeff;

				std::tie(bidx, coeff) = basis_.hamiltonianCoeff(s, n);

				if(bidx >= 0)
					m[bidx] += J2_*(1-sgn)*coeff;
			}
		}
		return m;
	}
};
#endif //CY_TIHAMXXZ_HPP
