#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <iostream>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include "edlib/Basis/Basis2D.hpp"
#include "edlib/Basis/FullBasis.hpp"
#include "edlib/Hamiltonians/TITFI2D.hpp"
#include "edlib/Op/NodeMV.hpp"

#include "edlib/EDP/LocalHamiltonian.hpp"
#include "edlib/EDP/ConstructSparseMat.hpp"

#include "utils.hpp"

using namespace edlib;
using namespace Eigen;

template<typename Basis>
Eigen::VectorXd translateY(Basis&& basis, const Eigen::VectorXd& r)
{
	Eigen::VectorXd res(r.size());
	for(int i = 0; i < r.size(); i++)
	{
		res(basis.rotateY(i,1)) = r(i);
	}
	return res;
}

template<typename Basis>
Eigen::VectorXd translateX(Basis&& basis, const Eigen::VectorXd& r)
{
	Eigen::VectorXd res(r.size());
	for(int i = 0; i < r.size(); i++)
	{
		res(basis.rotateX(i,1)) = r(i);
	}
	return res;
}

void TestBasisMatrix(const Eigen::MatrixXd& r)
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
void CheckBasis2DParity(Basis&& basis, const MatrixXd& r)
{
	const uint32_t kx = basis.getKx();
	const uint32_t ky = basis.getKy();
	const uint32_t p = basis.getP();

	const int expkx = (kx == 0)?1:-1;
	const int expky = (ky == 0)?1:-1;

	for(int i = 0; i < r.cols(); i++)
	{
		VectorXd c = r.col(i);
		REQUIRE( ( c - expkx*translateX(basis, c)).norm() < 1e-10);
		REQUIRE( ( c - expky*translateY(basis, c)).norm() < 1e-10);
		REQUIRE( ( c - p*flip(basis, c)).norm() < 1e-10);
	}
}

template<class Basis>
void CheckBasis2D(Basis&& basis, const MatrixXd& r)
{
	const uint32_t kx = basis.getKx();
	const uint32_t ky = basis.getKy();

	const int expkx = (kx == 0)?1:-1;
	const int expky = (ky == 0)?1:-1;

	for(int i = 0; i < r.cols(); i++)
	{
		VectorXd c = r.col(i);
		REQUIRE( ( c - expkx*translateX(basis, c)).norm() < 1e-10);
		REQUIRE( ( c - expky*translateY(basis, c)).norm() < 1e-10);
	}
}

Eigen::SparseMatrix<double> getSZZ()
{
	Eigen::SparseMatrix<double> res(4,4);
	res.insert(0,0) = 1.0;
	res.insert(1,1) = -1.0;
	res.insert(2,2) = -1.0;
	res.insert(3,3) = 1.0;
	res.makeCompressed();
	return res;
}

Eigen::SparseMatrix<double> getSX()
{
	Eigen::SparseMatrix<double> res(2,2);
	res.insert(0,1) = 1.0;
	res.insert(1,0) = 1.0;
	res.makeCompressed();
	return res;
}

TEST_CASE("Print", "[basis2d]")
{
	const int Lx = 3;
	const int Ly = 2;
	const int N = Lx*Ly;
	Basis2D<uint32_t> basis(Lx, Ly, 0, 0, false);
	
	Eigen::MatrixXd bmat = basisMatrix(basis);

	for(uint32_t n = 0; n < basis.getDim(); ++n)
	{
		std::cout << bmat.col(n).norm() << std::endl;
	}
}

TEST_CASE("Check Basis2D 4x4", "[basis2d]")
{
	const uint32_t Lx = 4;
	const uint32_t Ly = 4;

	SECTION("Not use U1")
	{
		for(const uint32_t kx: {0u, Lx/2})
		for(const uint32_t ky: {0u, Ly/2})
		{
			Basis2D<uint32_t> basis(Lx, Ly, kx, ky, false);
			printf("Test kx = %u, ky = %u, dim = %lu\n", kx, ky, basis.getDim());
			MatrixXd r = basisMatrix(basis);
			TestBasisMatrix(r);
			CheckBasis2D(basis, r);
		}
	}
	SECTION("Not use U1")
	{
		for(const uint32_t kx: {0u, Lx/2})
		for(const uint32_t ky: {0u, Ly/2})
		{
			Basis2D<uint32_t> basis(Lx, Ly, kx, ky, true);
			printf("Test kx = %u, ky = %u, dim = %lu\n", kx, ky, basis.getDim());
			MatrixXd r = basisMatrix(basis);
			TestBasisMatrix(r);
			CheckBasis2D(basis, r);
		}
	}
}

TEST_CASE("Compare enegies from the 2D TFI model", "[tfi-2d]") 
{
	const int Lx = 3;
	const int Ly = 2;
	const int N = Lx*Ly;
	Basis2D<uint32_t> basis2D(Lx, Ly, 0, 0, false);

	double h = 10.0;

	double gsEnergy1;
	{
		edp::LocalHamiltonian<double> lh(N, 2);
		for(uint32_t nx = 0; nx < Lx; ++nx)
		{
			for(uint32_t ny = 0; ny < Ly; ++ny)
			{
				lh.addTwoSiteTerm(
					std::make_pair(Lx*ny+nx, Lx*( (ny+1) % Ly)+nx), -getSZZ()
				);
				lh.addTwoSiteTerm(
					std::make_pair(Lx*ny+nx, Lx*ny+((nx+1) % Lx)), -getSZZ()
				);
				lh.addOneSiteTerm(Lx*ny+nx, -h*getSX());
			}
		}
		Eigen::SparseMatrix<double> m = edp::constructSparseMat<double>(1<<N, lh);

		Spectra::SparseSymMatProd<double> op(m);
		Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, decltype(op)> eigs(&op, 2, 6);
		eigs.init();
		eigs.compute(10000, 1e-12, Spectra::SMALLEST_ALGE);
		std::cout << eigs.eigenvalues() << std::endl;
		gsEnergy1 = eigs.eigenvalues()[0];
	}

	double gsEnergy2;
	{
		TITFI2D<uint32_t> ham(basis2D, 1.0, h);
		const int dim = basis2D.getDim();
		NodeMV mv(dim, 0, dim, ham);
		Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, NodeMV> eigs(&mv, 2, 4);
		eigs.init();
		eigs.compute(10000, 1e-12, Spectra::SMALLEST_ALGE);
		std::cout << eigs.eigenvalues() << std::endl;
		gsEnergy2 = eigs.eigenvalues()[0];
	}

	REQUIRE_THAT(gsEnergy1, Catch::Matchers::WithinAbs(gsEnergy2, 1e-8));
}
