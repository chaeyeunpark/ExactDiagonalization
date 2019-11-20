#include <iostream>
#include "Basis/TIBasis.hpp"
#include "Basis/TIBasisZ2.hpp"
#include "EDP/ConstructSparseMat.hpp"
#include "XXZ.hpp"
#include "Hamiltonians/TIXXZ.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

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

template<class Basis>
void CheckBasis(Basis&& basis)
{
	const int N = basis.getN();
	const int K = basis.getK();
	const int p = basis.getP();
	const int expk = (K==0)?1:-1;

	MatrixXd r = basis.basisMatrix().transpose();
	for(int i = 0; i < r.cols(); i++)
	{
		VectorXd c = r.col(i);
		REQUIRE( ( c - expk*translate(basis, c)).norm() < 1e-10);
		REQUIRE( ( c - p*flip(basis, c)).norm() < 1e-10);
		//std::cout << c.transpose() << std::endl;
	}

	TIXXZ<uint32_t> tiXXZ(basis, 1.0, 0.9);
	auto hamTI = edp::constructSparseMat<double>(basis.getDim(), [&tiXXZ](uint32_t col){return tiXXZ.getCol(col);});

	XXZ xxz(N, 1.0, 0.9);
	auto hamFull = edp::constructSparseMat<double>(1<<N, xxz);
	Eigen::MatrixXd mat = r.transpose()*hamFull*r;

	REQUIRE( (Eigen::MatrixXd(hamTI) - mat).norm() < 1e-10);
}

template<class Basis>
void PrintRpts(Basis&& basis)
{
	for(int i = 0; i < basis.getDim(); i++)
	{
		std::cout << basis.getNthRep(i) << std::endl;
	}
}
TEST_CASE("TIBasisZ2 test basis matrices", "[z2_basis_check]")
{
	const std::array<int, 2> ps{-1,1};
	
	for(int N = 4; N <= 12; N+=2)
	{
		for(int K: {0, N/2})
		{
			for(int p: ps)
			{
				printf("Testing N: %d, K: %d, p: %d\n", N, K, p);
				TIBasisZ2<uint32_t> basis(N,K,p);
				CheckBasis(basis);
			}
		}
	}
}
TEST_CASE("TIBasisZ2 test coeffs", "[z2_basis_check]")
{
	const std::array<int, 2> ps{-1,1};
	
	for(int N = 4; N <= 12; N+=2)
	{
		for(int K: {0, N/2})
		{
			for(int p: ps)
			{
				printf("Testing N: %d, K: %d, p: %d\n", N, K, p);
				TIBasisZ2<uint32_t> basis(N,K,p);

				auto mat = basis.basisMatrix();
				for(int s = 0; s < (1<<N); ++s)
				{
					int idx;
					double coeff;
					std::tie(idx, coeff) = basis.coeffAt(s);
					if(idx != -1)
					{
						REQUIRE(std::abs(coeff - mat(idx,s)) < 1e-7);
					}
				}
			}
		}
	}
}
