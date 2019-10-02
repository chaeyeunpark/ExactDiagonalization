#include <iostream>
#include <bitset>

#include <armadillo>

#include <SparseMat/ConstructMat.hpp>
#include <SparseMat/LocalSparseOperator.hpp>
#include <SparseMat/SpinMatrices.hpp>

#include "TIBasisLadder.hpp"
#include "TILadderYang.hpp"

class TestHelper
{
	ED::TIBasisLadder<uint32_t> basis_;

public:
	TestHelper(int n, int k, int parity, int px, int py)
		: basis_(n, k, parity, px, py)
	{
	}
	arma::cx_vec exhausitive(uint32_t sigma)
	{
		using namespace arma;
		constexpr std::complex<double> I(0., 1.);
		const int L = basis_.getN()/2;
		std::complex<double> expk = std::exp(2.0*I*M_PI*double(basis_.getK())/double(L));
		cx_vec res(1<<basis_.getN(), fill::zeros);

		for(int i = 0; i < 8; i++)
		{
			uint32_t s = basis_.reflect(sigma, i);
			int p = basis_.getP(i);
			for(int r = 0; r < basis_.getN()/2; r++)
			{
				res(basis_.rotl(s, 2*r)) += double(p)*std::pow(expk, r);
			}
		}
		return res;
	}
};

arma::mat FullHam(int N, double J, double g, double h)
{
	edp::LocalSparseOperator<double> lso(N, 2);
	const int L = N/2;

	arma::sp_mat xx(4,4);
	xx(0, 3) = 1.0;
	xx(1, 2) = 1.0;
	xx(2, 1) = 1.0;
	xx(3, 0) = 1.0;
	
	for(int i = 0; i < L; i++)
	{
		lso.addTwoSiteTerm(std::make_pair(2*i, 2*i+1), -g*xx);
	}
	for(int i = 0; i < N; i++)
	{
		lso.addOneSiteTerm(i, h*getX());
	}
	for(int i = 0; i < L; i++)
	{
		lso.addTwoSiteTerm(std::make_pair(2*i,(2*i+2)%N), -J*getSZZ());
		lso.addTwoSiteTerm(std::make_pair(2*i+1,(2*i+3)%N), -J*getSZZ());
	}
	return edp::constructMat<double>(1<<N, lso);
}

int main()
{
	using namespace ED;
	constexpr int N = 12;
	constexpr int K = 3;
	using UINT = uint32_t;
	TIBasisLadder<UINT> basis(N,K,1,1,1);
	TestHelper tt(N, K, 1,1,1);
	constexpr double J = 1.0;
	constexpr double g = 1.0;
	constexpr double h = 0.2;
	TILadderYang<UINT> ham(basis, J, g, h);

	arma::mat mat = edp::constructMat<double>(basis.getDim(), ham);

	std::cout << mat << std::endl;
	
	arma::cx_mat ex(1<<N, basis.getDim(), arma::fill::zeros);
	for(int i = 0; i < basis.getDim(); i++)
	{
		TIBasisLadder<uint32_t>::RepData t = basis.getRepData(i);
		ex.col(i) = arma::normalise(tt.exhausitive(basis.getNthRep(i)));
		std::cout << basis.getNthRep(i) << "\t" << t.rot << "\t";
		for(auto q: t.qs)
		{
			std::cout << q << "\t";
		}
		std::cout << std::endl;
	}

	arma::mat exReal = real(ex);

	arma::mat fullHam = FullHam(N, J, g, h);
	arma::mat redHam(basis.getDim(), basis.getDim());
	for(int i = 0; i < basis.getDim(); i++)
	{
		for(int j = 0; j < basis.getDim(); j++)
		{
			redHam(i,j) = dot(exReal.col(i), fullHam*exReal.col(j));
		}
	}

	std::cout << redHam << std::endl;

	assert(norm(redHam - mat) < 1e-6);

	return 0;
	
}
