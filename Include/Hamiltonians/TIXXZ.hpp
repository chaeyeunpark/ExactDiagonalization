#pragma once
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <map>
#include <boost/dynamic_bitset.hpp>

#include "Basis/Basis.hpp"

template<typename UINT>
class TIXXZ
{
private:
	const Basis<UINT>& basis_;
	double J_;
	double delta_;

public:
	TIXXZ(const Basis<UINT>& basis, double J, double delta)
		: basis_(basis), J_(J), delta_(delta)
	{
		
	}

	std::map<int, double> getCol(UINT n) const
	{
		unsigned int N = basis_.getN();

		UINT a = basis_.getNthRep(n);
		const boost::dynamic_bitset<> bs(N, a);

		std::map<int, double> m;
		for(unsigned int i = 0; i < N; i++)
		{
			//Next-nearest
			{
				unsigned int j = (i+1)%N;
				int sgn = (1-2*bs[i])*(1-2*bs[j]);

				m[n] += J_*delta_*sgn;
				
				UINT s = a;
				s ^= basis_.mask({i,j});

				int bidx;
				double coeff;

				std::tie(bidx, coeff) = basis_.hamiltonianCoeff(s, n);
				
				if(bidx >= 0)
				{
					m[bidx] += J_*(1.0-sgn)*coeff;
				}
			}
		}
		return m;
	}
};
