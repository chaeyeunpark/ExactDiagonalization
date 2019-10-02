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
class TILadderYang
{
private:
	Basis<UINT>& basis_;
	double J_;
	double g_;
	double h_;

public:
	TILadderYang(Basis<UINT>& basis, double J, double g, double h)
		: basis_(basis), J_(J), g_(g), h_(h)
	{
		
	}

	std::map<std::size_t,double> getCol(UINT aidx) const
	{
		const int N = basis_.getN();
		const int L = N/2;

		UINT sigma = basis_.getNthRep(aidx);
		const boost::dynamic_bitset<> bs(N, sigma);

		std::map<std::size_t, double> m;
		for(int n = 0; n < L; n++)
		{
			//-J*ZZ and hX for each alpha
			for(int alpha = 0; alpha < 2; alpha++)
			{
				int i = 2*n+alpha;
				int j = (2*n+2+alpha)%N;
				int sgn = (1-2*bs[i])*(1-2*bs[j]);

				m[aidx] += -J_*sgn;
				
				UINT s = sigma;
				UINT t = (1<<i);
				s ^= t;

				int bidx;
				double coeff;

				std::tie(bidx, coeff) = basis_.hamiltonianCoeff(s, aidx);
				
				if(bidx >= 0)
				{
					m[bidx] += h_*coeff;
				}
			}
			// XX for each ladder
			{
				UINT s = sigma;
				UINT t = (1 << 2*n) | (1 << (2*n + 1));
				s ^= t;
				int bidx;
				double coeff;
				std::tie(bidx, coeff) = basis_.hamiltonianCoeff(s, aidx);
				
				if(bidx >= 0)
				{
					m[bidx] += -g_*coeff;
				}
			}
		}
		return m;
	}
};
#endif //CY_TIHAMXXZ_HPP
