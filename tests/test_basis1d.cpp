#include <iostream>

#include "edlib/Basis/Basis1D.hpp"
#include "edlib/Basis/Basis1DZ2.hpp"
#include "edlib/Basis/ToOriginalBasis.hpp"
#include "edlib/EDP/ConstructSparseMat.hpp"
#include "edlib/Hamiltonians/TIXXZ.hpp"

#include "utils.hpp"
#include "XXZ.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace Eigen;
using namespace edlib;

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
void CheckBasis1DParity(Basis&& basis, const MatrixXd& r)
{
	const int N = basis.getN();
	const int k = basis.getK();
	const int p = basis.getP();

	const int expk = (k == 0)?1:-1;

	for(int i = 0; i < r.cols(); i++)
	{
		VectorXd c = r.col(i);
		REQUIRE( ( c - expk*translate(basis, c)).norm() < 1e-10);
		REQUIRE( ( c - p*flip(basis, c)).norm() < 1e-10);
	}
}

void TestBasisMatrix(const MatrixXd& r)
{
	using Catch::WithinAbs;
	for(int i = 0; i < r.cols(); i++)
	{
		auto c = r.col(i);
		REQUIRE_THAT( c.norm(), WithinAbs(1.0, 1e-6));
	}

	for(int i = 0; i < r.cols()-1; i++)
	{
		for(int j = i+1; j < r.cols(); j++)
		{
			auto c1 = r.col(i);
			auto c2 = r.col(j);
			REQUIRE_THAT( double(c1.transpose()*c2), WithinAbs(0.0, 1e-6));
		}
	}
}


template<class Basis>
void CheckBasis1D(Basis&& basis, const MatrixXd& r)
{
	const int N = basis.getN();
	const int k = basis.getK();
	const int expk = (k == 0)?1:-1;

	for(int i = 0; i < r.cols(); i++)
	{
		auto c = r.col(i);
		REQUIRE( ( c - expk*translate(basis, c)).norm() < 1e-10);
	}
}

template<class Basis>
void CheckBasisXXZ(Basis&& basis, const MatrixXd& r)
{
	TIXXZ<uint32_t> tiXXZ(basis, 1.0, 0.9);
	auto hamTI = edp::constructSparseMat<double>(basis.getDim(), 
			[&tiXXZ](uint32_t col){ return tiXXZ.getCol(col); });

	const int N = basis.getN();

	XXZ xxz(N, 1.0, 0.9);
	auto hamFull = edp::constructSparseMat<double>(1<<N, xxz);
	Eigen::MatrixXd mat = r.transpose()*hamFull*r;

	if( ! ((Eigen::MatrixXd(hamTI) - mat).norm() < 1e-10))
	{
		std::cout << mat << std::endl;
		std::cout << MatrixXd(hamTI) << std::endl;
	}
}

template<class Basis>
void PrintRpts(Basis&& basis)
{
	for(int i = 0; i < basis.getDim(); i++)
	{
		std::cout << basis.getNthRep(i) << std::endl;
	}
}

TEST_CASE("Basis1D works?", "[ti-basis]")
{
	int K = 0;
	for(int N = 20; N <= 28; N += 4)
	{
		Basis1D<uint32_t> basis(N, K, true);
		
		for(uint32_t i = 0; i < basis.getDim(); ++i)
		{
			uint32_t rep = basis.getNthRep(i);
			const auto [rep2, rot] = basis.findMinRots(rep);
			REQUIRE(rep2 == rep);
		}
	}
}

TEST_CASE("Basis1D test", "[basis]")
{
	SECTION("Not use U1")
	{
		for(int N = 4; N <= 12; N += 2)
		{
			for(int K: {0, N/2})
			{
				printf("Testing N: %d, K: %d\n", N, K);
				Basis1D<uint32_t> basis(N, K, false);
				MatrixXd r = basisMatrix(basis);
				TestBasisMatrix(r);
				CheckBasis1D(basis, r);
				CheckBasisXXZ(basis, r);
			}
		}
	}
	SECTION("Use U1")
	{
		for(int N = 4; N <= 12; N += 2)
		{
			for(int K: {0, N/2})
			{
				printf("Testing N: %d, K: %d\n", N, K);
				Basis1D<uint32_t> basis(N, K, true);
				MatrixXd r = basisMatrix(basis);
				TestBasisMatrix(r);
				CheckBasis1D(basis, r);
				CheckBasisXXZ(basis, r);
			}
		}
	}
}

TEST_CASE("Basis1DZ2 test", "[basis-z2]")
{
	const std::array<int, 2> ps{-1,1};

	SECTION("Not use U1")
	{
		for(int N = 4; N <= 12; N += 2)
		{
			for(int K: {0, N/2})
			{
				for(int p: ps)
				{
					printf("Testing N: %d, K: %d, p: %d\n", N, K, p);
					Basis1DZ2<uint32_t> basis(N, K, p, false);
					MatrixXd r = basisMatrix(basis);
					TestBasisMatrix(r);
					CheckBasis1DParity(basis, r);
					CheckBasisXXZ(basis, r);
				}
			}
		}
	}
	SECTION("Use U1")
	{
		for(int N = 4; N <= 12; N += 2)
		{
			for(int K: {0, N/2})
			{
				for(int p: ps)
				{
					printf("Testing N: %d, K: %d, p: %d\n", N, K, p);
					Basis1DZ2<uint32_t> basis(N, K, p, true);
					MatrixXd r = basisMatrix(basis);
					TestBasisMatrix(r);
					CheckBasis1D(basis, r);
					CheckBasisXXZ(basis, r);
				}
			}
		}
	}
}
