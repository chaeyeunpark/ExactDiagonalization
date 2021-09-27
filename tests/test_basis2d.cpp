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
