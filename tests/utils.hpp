#pragma once
#include <Eigen/Dense>

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
