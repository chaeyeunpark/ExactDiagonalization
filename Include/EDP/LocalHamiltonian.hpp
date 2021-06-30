#ifndef EDP_LOCALHAMILTONIAN_HPP
#define EDP_LOCALHAMILTONIAN_HPP
#include <Eigen/Sparse>
namespace edp
{

namespace internal
{
	inline uint32_t ipow(uint32_t base, uint32_t exp)
	{
		int result = 1u;
		while (exp)
		{
			if (exp & 1)
				result *= base;
			exp >>= 1;
			base *= base;
		}

		return result;
	}
}

template<typename T>
struct TwoSiteTerm
{
	std::pair<int, int> sites;
	Eigen::SparseMatrix<T> m;

	TwoSiteTerm(const std::pair<int,int>& p1, const Eigen::SparseMatrix<T>& p2)
		: sites(p1), m(p2)
	{
	}
	TwoSiteTerm(std::pair<int,int>&& p1, Eigen::SparseMatrix<T>&& p2)
		: sites(p1), m(p2)
	{
	}
};

template<typename T>
struct OneSiteTerm
{
	int site;
	Eigen::SparseMatrix<T> m;

	OneSiteTerm(int p1, const Eigen::SparseMatrix<T>& p2)
		: site(p1), m(p2)
	{
	}
	OneSiteTerm(int p1, Eigen::SparseMatrix<T>&& p2)
		: site(p1), m(p2)
	{
	}
};

template<typename T>
class LocalHamiltonian
{
private:
	int numSites_;
	int d_; //local Hibert dimension

	std::vector<TwoSiteTerm<T> > twoSiteTerms_;
	std::vector<OneSiteTerm<T> > oneSiteTerms_;

	int swapBaseD(int idx, int pos, int val) const
	{
		int b = internal::ipow(d_, pos);
		int upper = (idx / (b*d_))*d_ + val;
		return upper*b + (idx % b);
	}
public:
	LocalHamiltonian(int numSites, int d)
		: numSites_{numSites}, d_{d}
	{
	}
	LocalHamiltonian(const LocalHamiltonian&) = default;
	LocalHamiltonian(LocalHamiltonian&&) = default;

	LocalHamiltonian& operator=(const LocalHamiltonian&) = default;
	LocalHamiltonian& operator=(LocalHamiltonian&&) = default;

	void clearTerms()
	{
		std::vector<TwoSiteTerm<T> >().swap(twoSiteTerms_);
		std::vector<OneSiteTerm<T> >().swap(oneSiteTerms_);
	}
	int getNumSites() const
	{
		return numSites_;
	}

	std::map<uint32_t, T> getCol(uint32_t n) const;
	inline std::map<uint32_t, T> operator()(uint32_t n) const
	{
		return getCol(n);
	}

	void addTwoSiteTerm(const std::pair<int,int>& site, Eigen::SparseMatrix<T> m)
	{
		m.makeCompressed();
		twoSiteTerms_.emplace_back(site, std::move(m));
	}
	void addOneSiteTerm(int site, Eigen::SparseMatrix<T> m)
	{
		m.makeCompressed();
		oneSiteTerms_.emplace_back(site, std::move(m));
	}
};

}

template<typename T>
std::map<uint32_t, T> edp::LocalHamiltonian<T>::getCol(uint32_t n) const
{
	using internal::ipow;
	using Eigen::SparseMatrix;

	std::map<uint32_t, T> m;

	for(auto& twoSiteTerm : twoSiteTerms_)
	{
		auto a = (n/ipow(d_,twoSiteTerm.sites.first))%d_;
		auto b = (n/ipow(d_,twoSiteTerm.sites.second))%d_;
		//auto col = twoSiteTerm.m.col(a*d_ + b);
		auto col = b*d_ + a;
		for(typename SparseMatrix<T>::InnerIterator it(twoSiteTerm.m, col); it; ++it)
		{
			int r = it.row();
			int t = n;
			t = swapBaseD(t, twoSiteTerm.sites.first, r%d_);
			t = swapBaseD(t, twoSiteTerm.sites.second, r/d_);
			m[t] += it.value();
		}
	}
	for(auto& oneSiteTerm: oneSiteTerms_)
	{
		int a = (n/ipow(d_,oneSiteTerm.site))%d_;
		//auto col = oneSiteTerm.m.col(a);
		auto col = a;

		for(typename SparseMatrix<T>::InnerIterator it(oneSiteTerm.m, col); it; ++it)
		{
			int r = it.row();
			int t = n;
			t = swapBaseD(t, oneSiteTerm.site, r);
			m[t] += it.value();
		}
	}
	return m;
}


#endif//EDP_LOCALHAMILTONIAN_HPP
