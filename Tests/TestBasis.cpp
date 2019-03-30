#include <iostream>
#include "TIBasis.hpp"
#include "TIBasisZ2.hpp"
#include "EDP/ConstructSparseMat.hpp"
#include "XXZ.hpp"
#include "TIXXZ.hpp"
using namespace Eigen;
template<typename Basis>
VectorXd translate(Basis&& basis, const VectorXd& r)
{
	VectorXd res(r.size());
	for(int i = 0; i < r.size(); i++)
	{
		res(basis.rotl(i,1)) = r(i);
	}
	return res;
}

template<typename Basis>
VectorXd flip(Basis&& basis, const VectorXd& r)
{
	VectorXd res(r.size());
	for(int i = 0; i < r.size(); i++)
	{
		res(basis.flip(i)) = r(i);
	}
	return res;
}
int main(int argc, char* argv[])
{
	const int N = 8;
	const int K = 4;

	const int expk = (K==0)?1:-1;
	const int p = 1;

	TIBasisZ2<uint32_t> basis(N,K,p);

	for(int i = 0; i < basis.getDim(); i++)
	{
		std::cout << basis.getNthRep(i) << std::endl;
	}

	MatrixXd r = basis.basisMatrix().transpose();
	for(int i = 0; i < r.cols(); i++)
	{
		VectorXd c = r.col(i);
		std::cout << ( c - expk*translate(basis, c)).norm() << std::endl;
		std::cout << ( c - p*flip(basis, c)).norm() << std::endl;
		std::cout << c.transpose() << std::endl;
	}

	std::cout << "------------------------" << std::endl;

	TIXXZ<uint32_t> tiXXZ(basis, 1.0, 0.9);
	auto hamTI = edp::constructSparseMat<double>(basis.getDim(), [&tiXXZ](uint32_t col){return tiXXZ.getCol(col);});

	std::cout << Eigen::MatrixXd(hamTI) << std::endl;

	std::cout << "===========================" << std::endl;
	
	XXZ xxz(N, 1.0, 0.9);
	auto hamFull = edp::constructSparseMat<double>(1<<N, xxz);
	Eigen::MatrixXd mat = r.transpose()*hamFull*r;

	std::cout << mat << std::endl;

	std::cout << (Eigen::MatrixXd(hamTI) - mat).norm() << std::endl;

	
	/*
	MatrixXd hamCoeff(dim, 1<<N);
	const auto dim = basis.getDim();
	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < (1<<N); j++)
		{
			auto t = basis.hamiltonianCoeff(j, i);
			std::cout << t.second << "\t";
		}
		std::cout << std::endl;
	}
	*/
	return 0;
}
