#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <Spectra/SymEigsSolver.h>
#include <fstream>
#include <ios>

//#include "TIBasis.hpp"
#include "TIBasisZ2.hpp"
#include "TITFIsing.hpp"
#include "NodeMV.hpp"

int main(int argc, char* argv[])
{
	constexpr int N = 28;
	using UINT = uint32_t;

	std::cout << "#N: " << N << std::endl;

	std::vector<double> hs;
	for(int i = 0; i <= 20; i++)
	{
		hs.emplace_back(i*0.1);
	}
	for(auto h : hs)
	{
		Eigen::VectorXd ev[2];
		{
			TIBasisZ2<UINT> basis(N, 0, 1);
			TITFIsing<UINT> ham(basis, 1.0, h);
			const int dim = basis.getDim();

			NodeMV mv(dim, 0, dim, ham);

			Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, NodeMV> eigs(&mv, 2, 6);
			eigs.init();
			eigs.compute(10000, 1e-12, Spectra::SMALLEST_ALGE);
			if(eigs.info() != Spectra::SUCCESSFUL)
				return 1;
			ev[0] = eigs.eigenvalues();
		}
		{
			TIBasisZ2<UINT> basis(N, 0, -1);
			TITFIsing<UINT> ham(basis, 1.0, h);
			const int dim = basis.getDim();

			NodeMV mv(dim, 0, dim, ham);

			Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, NodeMV> eigs(&mv, 2, 6);
			eigs.init();
			eigs.compute(10000, 1e-12, Spectra::SMALLEST_ALGE);
			if(eigs.info() != Spectra::SUCCESSFUL)
				return 1;
			ev[1] = eigs.eigenvalues();
		}

		printf("%f\t%.10f\t%.10f\t%.10f\t%.10f\n",h,ev[0](0), ev[0](1), ev[1](0), ev[1](1));
	}

	return 0;
}
