#pragma once
#include <Eigen/Dense>

inline void make_first_positive(Eigen::VectorXi& v)
{
	for(uint32_t n = 0; n < v.size(); ++ n) 
	{
		if(v(n) == 0)
			continue;
		if(v(n) < 0)
		{
			v *= -1;
			return ;
		}
		else
			return ;
	}
}


struct hash_vector
{
	size_t operator()(const Eigen::VectorXi& v) const
	{
		std::size_t seed = v.size();
		for(uint32_t n = 0; n < v.size(); ++ n) {
			seed ^= abs(v(n)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

int powi(int base, unsigned int exp)
{
    int res = 1;
    while (exp) {
        if (exp & 1)
            res *= base;
        exp >>= 1;
        base *= base;
    }
    return res;
}

template<typename UINT, template<typename> class Basis>
Eigen::MatrixXd basisMatrix(const Basis<UINT>& basis)
{
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(1u << basis.getN(), basis.getDim());
	for(unsigned int n = 0; n < basis.getDim(); ++n)
	{
		auto bvec = basis.basisVec(n);
		for(const auto p: bvec)
		{
			res(p.first, n) = p.second;
		}
	}
	return res;
}

template<typename Basis>
Eigen::VectorXd flip(Basis&& basis, const Eigen::VectorXd& r)
{
	Eigen::VectorXd res(r.size());
	for(int i = 0; i < r.size(); i++)
	{
		res(basis.flip(i)) = r(i);
	}
	return res;
}


