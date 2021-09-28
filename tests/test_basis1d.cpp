#include <iostream>
#include <unordered_set>

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

constexpr uint32_t MAX_N = 14;

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

template<class Basis>
uint64_t countDim1D(Basis&& basis, bool useU1)
{
	const int N = basis.getN();
	const int k = basis.getK();
	const int expk = (k == 0)?1:-1;

	std::unordered_set<VectorXi, hash_vector> basisVecs;

	for(uint32_t n = 0; n < (1u << N); ++ n)
	{
		if(useU1 && (__builtin_popcountll(n) != N/2))
			continue;
		Eigen::VectorXi v = Eigen::VectorXi::Zero(1ul << N);
		for(uint32_t k = 0; k < N; k++)
		{
			v(basis.rotl(n, k)) += powi(expk, k);
		}
		if (v.cwiseAbs().sum() != 0)
		{
			make_first_positive(v);
			basisVecs.emplace(std::move(v));
		}
	}

	return basisVecs.size();
}

template<class Basis>
uint64_t countDim1DZ2(Basis&& basis, bool useU1)
{
	const int N = basis.getN();
	const int k = basis.getK();
	const int p = basis.getP();
	const int expk = (k == 0)?1:-1;

	std::unordered_set<VectorXi, hash_vector> basisVecs;

	for(uint32_t n = 0; n < (1u << N); ++ n)
	{
		if(useU1 && (__builtin_popcountll(n) != N/2))
			continue;
		Eigen::VectorXi v = Eigen::VectorXi::Zero(1ul << N);
		for(uint32_t k = 0; k < N; k++)
		{
			v(basis.rotl(n, k)) += powi(expk, k);
			v(basis.rotl(basis.flip(n), k)) += powi(expk, k)*p;
		}
		if (v.cwiseAbs().sum() != 0)
		{
			make_first_positive(v);
			basisVecs.emplace(std::move(v));
		}
	}

	return basisVecs.size();
}

void TestBasisMatrix(const Eigen::MatrixXd& r)
{
	using Catch::WithinAbs;
	Eigen::MatrixXd id = Eigen::MatrixXd::Identity(r.cols(), r.cols());
	REQUIRE_THAT((r.transpose()*r - id).cwiseAbs().maxCoeff(),
			WithinAbs(0.0, 1e-8));
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

	REQUIRE((Eigen::MatrixXd(hamTI) - mat).norm() < 1e-10);
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

TEST_CASE("Basis1D test", "[basis1d]")
{
	SECTION("Not use U1")
	{
		for(int N = 4; N <= MAX_N; N += 2)
		{
			for(int K: {0, N/2})
			{
				printf("Testing N: %d, K: %d\n", N, K);
				Basis1D<uint32_t> basis(N, K, false);

				REQUIRE(basis.getDim() == countDim1D(basis, false));

				MatrixXd r = basisMatrix(basis);
				TestBasisMatrix(r);
				CheckBasis1D(basis, r);
				CheckBasisXXZ(basis, r);
			}
		}
	}
	SECTION("Use U1")
	{
		for(int N = 4; N <= MAX_N; N += 2)
		{
			for(int K: {0, N/2})
			{
				printf("Testing N: %d, K: %d\n", N, K);
				Basis1D<uint32_t> basis(N, K, true);

				REQUIRE(basis.getDim() == countDim1D(basis, true));

				MatrixXd r = basisMatrix(basis);
				TestBasisMatrix(r);
				CheckBasis1D(basis, r);
				CheckBasisXXZ(basis, r);
			}
		}
	}
}

TEST_CASE("Basis1DZ2 test", "[basis1dz2]")
{
	const std::array<int, 2> ps{-1,1};
	
	SECTION("Not use U1")
	{
		for(int N = 4; N <= MAX_N; N += 2)
		{
			for(int K: {0, N/2})
			{
				for(int p: ps)
				{
					printf("Testing N: %d, K: %d, p: %d\n", N, K, p);
					Basis1DZ2<uint32_t> basis(N, K, p, false);

					REQUIRE(basis.getDim() == countDim1DZ2(basis, false));

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
		for(int N = 6; N <= MAX_N; N += 2)
		{
			for(int K: {0, N/2})
			{
				for(int p: ps)
				{
					printf("Testing N: %d, K: %d, p: %d\n", N, K, p);
					Basis1DZ2<uint32_t> basis(N, K, p, true);

					REQUIRE(basis.getDim() == countDim1DZ2(basis, true));

					MatrixXd r = basisMatrix(basis);
					TestBasisMatrix(r);
					CheckBasis1D(basis, r);
					CheckBasisXXZ(basis, r);
				}
			}
		}
	}
}
