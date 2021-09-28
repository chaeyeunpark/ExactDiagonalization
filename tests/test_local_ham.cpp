#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <unsupported/Eigen/KroneckerProduct>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include <iostream>
#include <cassert>
#include <random>
#include <algorithm>

#include "edlib/EDP/LocalHamiltonian.hpp"
#include "edlib/EDP/ConstructSparseMat.hpp"

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
	const uint32_t dim = (1<<N);

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
	const uint32_t dim = (1<<N);

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

class GraphGenerator
{
private:
	int numVertices_;
	std::vector<std::pair<int,int> > allEdges_;
	std::vector<int> allVertices_;

public:
	GraphGenerator(int numVertices)
		: numVertices_{numVertices}
	{
		for(int i = 0; i < numVertices-1; ++i)
		{
			for(int j = i+1; j < numVertices; ++j)
			{
				allEdges_.emplace_back(i,j);
			}
		}
		for(int i = 0; i < numVertices_; ++i)
		{
			allVertices_.emplace_back(i);
		}
	}

	std::vector<std::pair<int,int> > createRandomGraph(int numEdges)
	{
		std::vector<std::pair<int,int> > edges = allEdges_;
		std::random_shuffle(edges.begin(), edges.end());
		edges.resize(numEdges);
		return edges;
	}

	std::vector<int> createRandomVertexSet(int n)
	{
		std::vector<int> v = allVertices_;
		std::random_shuffle(v.begin(), v.end());
		v.resize(n);
		return v;
	}
};

constexpr int N = 10;
TEST_CASE("Test single qubit operators", "[LocalHamSingle]") {
	std::random_device rd;
	std::default_random_engine re{rd()};
	std::uniform_int_distribution<> uid(0, N-1);
	std::uniform_real_distribution<> urd;

	using cx_double = std::complex<double>;

	SECTION("Test contructing pauli X") {
		edp::LocalHamiltonian<double> lh(N,2);
		for(int n = 0; n < 100; ++n)
		{
			lh.clearTerms();
			int idx = uid(re);
			double val = urd(re);
			lh.addOneSiteTerm(idx, val*getSX());
			auto mat1 = Eigen::MatrixXd(edp::constructSparseMat<double>(1<<N, lh));
			auto mat2 = singleQubitOp<double>(N, idx, val*getSX()) ;
			REQUIRE((mat1-mat2).squaredNorm() < 1e-8);
		}
	}
	SECTION("Test contructing pauli Z") {
		edp::LocalHamiltonian<double> lh(N,2);
		for(int n = 0; n < 100; ++n)
		{
			lh.clearTerms();
			int idx = uid(re);
			double val = urd(re);
			lh.addOneSiteTerm(idx, val*getSZ());
			auto mat1 = Eigen::MatrixXd(edp::constructSparseMat<double>(1<<N, lh));
			auto mat2 = singleQubitOp<double>(N, idx, val*getSZ()) ;
			REQUIRE((mat1-mat2).squaredNorm() < 1e-8);
		}
	}
	SECTION("Test contructing pauli Y") {
		edp::LocalHamiltonian<cx_double> lh(N,2);
		for(int n = 0; n < 100; ++n)
		{
			lh.clearTerms();
			int idx = uid(re);
			double val = urd(re);
			lh.addOneSiteTerm(idx, val*getSY());
			auto mat1 = Eigen::MatrixXcd(edp::constructSparseMat<cx_double>(1<<N, lh));
			auto mat2 = singleQubitOp<cx_double>(N, idx, val*getSY()) ;
			REQUIRE((mat1-mat2).squaredNorm() < 1e-8);
		}
	}
}

TEST_CASE("Test 2-local Hamiltonians", "[LocalHam2loc]") {
	std::random_device rd;
	std::default_random_engine re{rd()};
	std::uniform_int_distribution<> uid(0, N-1);
	std::uniform_int_distribution<> numEdgesRd(1, N*(N-1)/2);
	std::uniform_real_distribution<> urd;

	using cx_double = std::complex<double>;

	GraphGenerator ggen{N};

	SECTION("Random XYZ Hamiltonian") {
		edp::LocalHamiltonian<double> lh(N,2);

		for(int n = 0; n < 10; ++n)
		{
			lh.clearTerms();
			int numEdges = numEdgesRd(re);
			auto edges = ggen.createRandomGraph(numEdges);

			auto matEx = Eigen::MatrixXd(1<<N, 1<<N);
			matEx.setZero();

			for(const auto& edge: edges)
			{
				double v1 = urd(re);
				double v2 = urd(re);
				double v3 = urd(re);
				lh.addTwoSiteTerm(edge, v1*getSXX());
				lh.addTwoSiteTerm(edge, v2*getSYY());
				lh.addTwoSiteTerm(edge, v3*getSZZ());

				matEx += v1*twoQubitOp<double>(N, edge.first, edge.second, getSX(), getSX());
				auto matYY = twoQubitOp<cx_double>(N, edge.first, edge.second, getSY(), getSY());
				matEx += v2*matYY.real();
				matEx += v3*twoQubitOp<double>(N, edge.first, edge.second, getSZ(), getSZ());
			}
			auto mat = Eigen::MatrixXd(edp::constructSparseMat<double>(1<<N, lh));
			REQUIRE((mat-matEx).squaredNorm() < 1e-8);
		}
	}
	SECTION("Random Field Ising Hamiltonian") {
		edp::LocalHamiltonian<double> lh(N,2);

		for(int n = 0; n < 20; ++n)
		{
			lh.clearTerms();
			int numEdges = numEdgesRd(re);
			auto edges = ggen.createRandomGraph(numEdges);

			auto matEx = Eigen::MatrixXd(1<<N, 1<<N);
			matEx.setZero();

			for(const auto& edge: edges)
			{
				double v = urd(re);
				lh.addTwoSiteTerm(edge, v*getSZZ());

				matEx += v*twoQubitOp<double>(N, edge.first, edge.second, getSZ(), getSZ());
			}
			
			auto numFields = uid(re);
			auto vs = ggen.createRandomVertexSet(numFields);

			for(auto v: vs)
			{
				double h1 = urd(re);
				double h3 = urd(re);

				lh.addOneSiteTerm(v, h1*getSX());
				lh.addOneSiteTerm(v, h3*getSZ());

				matEx += h1*singleQubitOp<double>(N, v, getSX());
				matEx += h3*singleQubitOp<double>(N, v, getSZ());
			}

			auto mat = Eigen::MatrixXd(edp::constructSparseMat<double>(1<<N, lh));
			REQUIRE((mat-matEx).squaredNorm() < 1e-8);
		}
	}
}
