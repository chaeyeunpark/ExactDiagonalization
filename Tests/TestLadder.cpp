#include <iostream>
#include <bitset>
#include "TIBasisLadder.hpp"

#include <Eigen/Dense>

class TestHelper
{
	ED::TIBasisLadder<uint32_t> basis_;

public:
	TestHelper(int n, int k, int parity, int px, int py)
		: basis_(n, k, parity, px, py)
	{
	}
	Eigen::VectorXcd exhausitive(uint32_t sigma)
	{
		using namespace Eigen;
		constexpr std::complex<double> I(0., 1.);
		const int L = basis_.getN()/2;
		std::complex<double> expk = std::exp(2.0*I*M_PI*double(basis_.getK())/double(L));
		VectorXcd res = VectorXcd::Zero(1<<basis_.getN());

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

int main()
{
	using namespace ED;
	constexpr int N = 12;
	constexpr int K = 0;
	TIBasisLadder<uint32_t> basis(N,K,1,1,-1);
	TestHelper tt(N,K,1,1,-1);

	std::cout << basis.getDim() << std::endl;

	Eigen::MatrixXcd ex(1<<N, 1<<N);
	ex.setZero();
	for(int i = 0; i < basis.getDim(); i++)
	{
		TIBasisLadder<uint32_t>::RepData t = basis.getRepData(i);
		std::cout << basis.getNthRep(i) << "\t" << t.rot << "\t";
		for(auto q : t.qs)
			std::cout << q << "\t";
		std::cout << std::endl;
		ex.col(basis.getNthRep(i)) = tt.exhausitive(basis.getNthRep(i)).normalized();
	}
	Eigen::MatrixXcd ex1(1<<N, 1<<N);
	for(uint32_t s = 0; s < (1u << N); ++s)
	{
		auto t = basis.coeff(s);
		if(t.first < 0)
			continue;
		ex1(s, basis.getNthRep(t.first)) = t.second;
	}

	assert( (ex1 - ex).norm() < 1e-6);
}
