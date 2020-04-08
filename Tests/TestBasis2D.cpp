#include <iostream>
#include "Basis/TIBasis2D.hpp"
#include "Basis/FullBasis.hpp"
#include "Hamiltonians/TITFI2D.hpp"
#include "NodeMV.hpp"

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include <EDP/LocalHamiltonian.hpp>
#include <EDP/ConstructSparseMat.hpp>

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

int main()
{
	const int Lx = 3;
	const int Ly = 2;
	const int N = Lx*Ly;
	TIBasis2D<uint32_t> basis2D(Lx, Ly, 0 ,0, false);

	for(int i = 0; i < basis2D.getDim(); i++)
	{
		std::cout << basis2D.getNthRep(i) << std::endl;
		for(const auto p: basis2D.basisVec(i))
		{
			std::cout << p.first << "\t" 
				<< p.second << std::endl;
		}
	}

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
				lh.addOneSiteTerm(Lx*ny+nx, -3.0*getSX());
			}
		}
		Eigen::SparseMatrix<double> m = edp::constructSparseMat<double>(1<<N, lh);
		/*
		Eigen::MatrixXd hamFull = m;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
		es.compute(hamFull);

		std::cout << es.eigenvalues() << std::endl;
		*/

		Spectra::SparseSymMatProd<double> op(m);
		Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, decltype(op)> eigs(&op, 2, 6);
		eigs.init();
		eigs.compute(10000, 1e-12, Spectra::SMALLEST_ALGE);
		std::cout << eigs.eigenvalues() << std::endl;
		gsEnergy1 = eigs.eigenvalues()[0];
	}

	double gsEnergy2;
	{
		TITFI2D<uint32_t> ham(basis2D, 1.0, 3.0);
		const int dim = basis2D.getDim();
		NodeMV mv(dim, 0, dim, ham);
		Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, NodeMV> eigs(&mv, 2, 4);
		eigs.init();
		eigs.compute(10000, 1e-12, Spectra::SMALLEST_ALGE);
		std::cout << eigs.eigenvalues() << std::endl;
		gsEnergy2 = eigs.eigenvalues()[0];
	}

	std::cout << gsEnergy1 << "\t" << gsEnergy2 << std::endl;


	return 0;
}
