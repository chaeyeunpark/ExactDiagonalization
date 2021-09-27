#pragma once
#include <Eigen/Eigen>

class XXZ
{
private:
	int n_;
	double J_;
	double Delta_;

public:
	XXZ(int n, double J, double Delta)
		: n_(n), J_(J), Delta_(Delta)
	{
	}
	
	template<class State>
	typename State::T operator()(const State& smp) const
	{
		typename State::T s = 0.0;
		//Nearest-neighbor
		for(int i = 0; i < n_; i++)
		{
			int zz = smp.sigmaAt(i)*smp.sigmaAt((i+1)%n_);
			s += J_*Delta_*zz; //zz
			s += J_*(1-zz)*smp.ratio(i, (i+1)%n_); //xx+yy
		}
		return s;
	}

	std::vector<std::array<int,2> > offDiagonals(const Eigen::VectorXi& s) const
	{
		std::vector< std::array<int, 2> > res;
		for(int i = 0; i < n_; i++)
		{
			if( s(i) != s((i+1)%n_) )
				res.push_back(std::array<int,2>{ i,(i+1)%n_ });
		}
		return res;
	}

	std::map<uint32_t, double> operator()(uint32_t col) const
	{
		std::map<uint32_t, double> m;
		for(int i = 0; i < n_; i++)
		{
			int b1 = (col >> i) & 1;
			int b2 = (col >> ((i+1)%n_)) & 1;
			int zz = (1-2*b1)*(1-2*b2);
			long long int x = (1 << i) | (1 << ((i+1)%(n_)));
			m[col] += J_*Delta_*zz;
			m[col ^ x] += J_*(1 - zz);
		}
		return m;
	}
};

