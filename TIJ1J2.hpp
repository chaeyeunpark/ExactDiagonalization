#ifndef CY_TIXXZ_HPP
#define CY_TIXXZ_HPP
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <map>
#include <boost/dynamic_bitset.hpp>
//#include "BitOperations.h"
#include "Basis.hpp"

template<typename UINT>
class TIXXZ
{
private:
	Basis<UINT>& basis_;
	double J1_;
	double J2_;

public:
	TIXXZ(Basis<UINT>& basis, double J1, double J2)
		: basis_(basis), J1_(J1), J2_(J2)
	{
		
	}

	std::map<std::size_t,double> getCol(UINT n) const
	{
		int N = basis_.getN();

		UINT a = basis_.getNthRep(n);
		const boost::dynamic_bitset<> bs(N, a);

		std::map<std::size_t, double> m;
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
				{
					m[bidx] += J1_*(1.0-sgn)*coeff;
				}
			}
			//Next-nearest
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
				{
					m[bidx] += J2_*(1.0-sgn)*coeff;
				}
			}
		}
			
		return m;
	}
};
#endif //CY_TIHAMXXZ_HPP
