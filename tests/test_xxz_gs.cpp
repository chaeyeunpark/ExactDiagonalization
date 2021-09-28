#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues> 

#include <unsupported/Eigen/KroneckerProduct>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>


#include <iostream>
#include <cassert>
#include <random>
#include <algorithm>

#include "edlib/Basis/Basis1DZ2.hpp"
#include "edlib/Basis/ToOriginalBasis.hpp"
#include "edlib/Hamiltonians/TIXXZ.hpp"
#include "edlib/Op/NodeMV.hpp"

#include "edlib/EDP/LocalHamiltonian.hpp"
#include "edlib/EDP/ConstructSparseMat.hpp"

using namespace edlib;

Eigen::SparseMatrix<double> getSX()
{
	Eigen::SparseMatrix<double> res(2,2);
	res.insert(0,1) = 1.0;
	res.insert(1,0) = 1.0;
	res.makeCompressed();
	return res;
}

Eigen::SparseMatrix<std::complex<double> > getSY()
{
	Eigen::SparseMatrix<std::complex<double> > res(2,2);
	constexpr std::complex<double> I(0., 1.);
	res.insert(0,1) = -I;
	res.insert(1,0) = I;
	res.makeCompressed();
	return res;
}

Eigen::SparseMatrix<double> getSZ()
{
	Eigen::SparseMatrix<double> res(2,2);
	res.insert(0,0) = 1.0;
	res.insert(1,1) = -1.0;
	res.makeCompressed();
	return res;
}

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
twoQubitOp(int N, int pos1, int pos2,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& v1,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& v2)
{
	using namespace Eigen;
	const uint32_t dim = (1u << N);

	assert(pos1 < pos2);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> res(1,1);
	res(0,0) = 1.0;

	for(int i = 0; i < pos1; i++)
	{
		res = Eigen::kroneckerProduct(MatrixXd::Identity(2,2), res).eval();
	}

	res = Eigen::kroneckerProduct(v1, res).eval();
	for(int i = pos1+1; i < pos2; i++)
	{
		res = Eigen::kroneckerProduct(MatrixXd::Identity(2,2), res).eval();
	}
	res = Eigen::kroneckerProduct(v2, res).eval();
	for(int i = pos2+1; i < N; i++)
	{
		res = Eigen::kroneckerProduct(MatrixXd::Identity(2,2), res).eval();
	}
	return res;
}

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
singleQubitOp(int N, int pos, 
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& v)
{
	using namespace Eigen;
	const uint32_t dim = (1u << N);

	using MatrixT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ;
	MatrixT res(1,1);
	res(0,0) = 1.0;

	for(int i = 0; i < pos; i++)
	{
		res = Eigen::kroneckerProduct(MatrixT::Identity(2,2), res).eval();
	}

	res = Eigen::kroneckerProduct(v, res).eval();
	for(int i = pos+1; i < N; i++)
	{
		res = Eigen::kroneckerProduct(MatrixT::Identity(2,2), res).eval();
	}
	return res;
}

Eigen::SparseMatrix<double> getSXXYY()
{
	Eigen::SparseMatrix<double> res(4,4);
	res.insert(1,2) = 2.0;
	res.insert(2,1) = 2.0;
	res.makeCompressed();
	return res;
}

Eigen::SparseMatrix<double> getSXX()
{
	Eigen::SparseMatrix<double> res(4,4);
	res.insert(0,3) = 1.0;
	res.insert(1,2) = 1.0;
	res.insert(2,1) = 1.0;
	res.insert(3,0) = 1.0;
	res.makeCompressed();
	return res;
}

Eigen::SparseMatrix<double> getSYY()
{
	Eigen::SparseMatrix<double> res(4,4);
	res.insert(0,3) = -1.0;
	res.insert(1,2) = 1.0;
	res.insert(2,1) = 1.0;
	res.insert(3,0) = -1.0;
	res.makeCompressed();
	return res;
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

template<int N>
class CompareXXZ
{
private:
	double delta_;
	Basis1DZ2<uint32_t> basis_;

	Eigen::MatrixXd hamFull_;
	double gsEnergy_;
	Eigen::VectorXd gsVec_;
public:

	constexpr static int k = (N/2)*((N/2)%2);
	constexpr static int parity = 1-2*((N/2)%2);

	CompareXXZ(double delta)
		: delta_{delta}, basis_{N, k, parity, true}
	{
		using namespace Eigen;
		static_assert(N % 2 == 0, "N must be even");

		edp::LocalHamiltonian<double> lh(N,2);
		for(int i = 0; i < N; i++)
		{
			lh.addTwoSiteTerm({i, (i+1) % N}, getSXXYY() + delta_*getSZZ());
		}
		hamFull_ = MatrixXd(edp::constructSparseMat<double>(1<<N, lh));
		SelfAdjointEigenSolver<MatrixXd> es;
		es.compute(hamFull_);

		gsEnergy_ = es.eigenvalues()[0];
		gsVec_ = es.eigenvectors().col(0);

	}

	void Test()
	{
		using namespace Eigen;
		TIXXZ<uint32_t> ham(basis_, 1.0, delta_);
		const int dim = basis_.getDim();

		NodeMV mv(dim, 0, dim, ham);

		Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, NodeMV> eigs(&mv, 2, 6);
		eigs.init();
		eigs.compute(10000, 1e-12, Spectra::SMALLEST_ALGE);
		if(eigs.info() != Spectra::SUCCESSFUL)
			REQUIRE(false);
		double gsEnergy1 = eigs.eigenvalues()[0];

		REQUIRE(std::abs(gsEnergy_ - gsEnergy1) < 1e-4);
		VectorXd subspaceGs = eigs.eigenvectors().col(0);
		VectorXd gsVec1;
		{
			auto v = toOriginalVector(basis_, subspaceGs.data());
			gsVec1 = Map<VectorXd>(v.data(), 1<<N);
		}

		double gsEnergy2 = double(gsVec1.transpose()*hamFull_*gsVec1)/double(gsVec1.transpose()*gsVec1);

		REQUIRE(std::abs(gsEnergy1 - gsEnergy2) < 1e-6);
		REQUIRE(std::abs(std::abs(gsVec_.transpose()*gsVec1) - 1.) < 1e-6);
	}
};

TEST_CASE("Compare GS of XXZ using LocalHamiltonian and TIBasis", "[XXZGS]") 
{
	SECTION("TIBasis Z2 XXZ N=8") {
		CompareXXZ<8> test(1.0);
		test.Test();
	}
	SECTION("TIBasis Z2 XXZ N=10") {
		CompareXXZ<10> test(1.0);
		test.Test();
	}
	SECTION("TIBasis Z2 XXZ N=12") {
		CompareXXZ<12> test(1.0);
		test.Test();
	}
}
