#pragma once
#include <cassert>

#include "Basis.hpp"

template<typename UINT>
class Basis2D
	: public Basis<UINT>
{
protected:
	const uint32_t Lx_;
	const uint32_t Ly_;
	const uint32_t kx_;
	const uint32_t ky_;
	const UINT px_ = (~UINT(0)) >> (sizeof(UINT)*8 - Lx_);

	/**
	 * In two dimension, there are several different rotations that gives 
	 * the same representation
	 * */
	std::tuple<UINT, uint32_t, uint32_t> getMinRots(UINT sigma) const
	{
		UINT rep = sigma;
		uint32_t rotX = Lx_;
		uint32_t rotY = Ly_;
		for(uint32_t rx = 1; rx <= Lx_; rx++)
		{
			UINT srx = rotateX(sigma, rx);
			for(uint32_t ry = 1; ry <= Ly_; ry++)
			{
				UINT sr = rotateY(srx, ry);
				if(sr < rep)
				{
					rep = sr;
					rotX = rx;
					rotY = ry;
				}
			}
		}
		return std::make_tuple(rep, rotX, rotY);
	}

	uint32_t checkState(UINT s) const
	{
		uint32_t cnt = 0;
		UINT sr = s;

		for(uint32_t rx = 1; rx <= Lx_; rx++)
		{
			for(uint32_t ry = 1; ry <= Ly_; ry++)
			{
				sr = rotateX(s, rx);
				sr = rotateY(sr, ry);
				if(sr < s)
				{
					return 0; //s is not a representative
				}
				else if(sr == s)
				{
					if(phase(rx,ry) == -1) //not allowed
					{
						return 0;
					}
					else
					{
						++cnt;
					}
				}
			}
		}
		return cnt;
	}

	int phase(int32_t rotX, int32_t rotY) const
	{
		//return exp(2\Pi I(kx_*rotX/Lx_ + ky_*rotY/Ly_))
		int sgn = 1;
		if((kx_*rotX) % Lx_ != 0)
		{
			sgn *= -1;
		}
		if((ky_*rotY) % Ly_ != 0)
		{
			sgn *= -1;
		}
		return sgn;
	}

public:
	Basis2D(uint32_t Lx, uint32_t Ly, uint32_t kx, uint32_t ky)
		: Basis<UINT>(Lx*Ly), Lx_{Lx}, Ly_{Ly}, kx_{kx}, ky_{ky}
	{
	}

	UINT rotateY(UINT sigma, int32_t r) const
	{
		return this->rotl(sigma, Lx_*r);
	}

	UINT rotateX(UINT sigma, int32_t r) const
	{
		assert((r >= 0) && (r <= Lx_));
		const auto Lx = Lx_;
		const auto Ly = Ly_;
		const auto px = px_;
		for(uint32_t ny = 0; ny < Ly; ny++)
		{
			UINT t = (sigma >> (ny*Lx)) & px;
			t = ((t << r) | (t >> (Lx-r))) & px;
			sigma &= ~(px << ny*Lx);
			sigma |= (t << ny*Lx);
		}
		return sigma;
	}

	inline uint32_t getLx() const { return Lx_; }
	inline uint32_t getLy() const { return Ly_; }
	inline uint32_t getKx() const { return kx_; }
	inline uint32_t getKy() const { return ky_; }


	inline uint32_t toIdx(uint32_t nx, uint32_t ny) const 
	{
		return ny*Lx_ + nx;
	}

};
