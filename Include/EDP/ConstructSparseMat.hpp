#ifndef EDP_CONSTRUCTSPARSEMAT_HPP
#define EDP_CONSTRUCTSPARSEMAT_HPP
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace edp
{
	template<typename T, class ColFunc>
	auto constructMat(uint64_t dim, ColFunc&& colFunc)
	-> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> res(dim, dim);
		res.setZero();

		for(uint64_t i = 0; i < dim; i++)
		{
			auto m = colFunc(i);
			for(auto& elt: m)
			{
				res(elt.first, i) = elt.second;
			}
		}
		return res;
	}

	template<typename T, class ColFunc>
	auto constructSparseMat(uint64_t dim, ColFunc&& colFunc)
	-> Eigen::SparseMatrix<T> 
	{
		using TripletT = Eigen::Triplet<T>;
		std::vector<TripletT> tripletList;
		tripletList.reserve(3*dim);
		for(uint64_t col = 0; col < dim; ++col)
		{
			auto m = colFunc(col);
			for(const auto& v: m)
			{
				tripletList.emplace_back(v.first, col, v.second);
			}
		}

		Eigen::SparseMatrix<T> res(dim, dim);
		res.setFromTriplets(tripletList.begin(), tripletList.end());
		return res;
	}

	//basis must be sorted
	template<typename T, typename ColFunc>
	auto constructSubspaceMat(ColFunc&& t, const std::vector<uint32_t>& basis) 
	-> Eigen::SparseMatrix<T> 
	{
		const int n = basis.size();

		using TripletT = Eigen::Triplet<T>;
		std::vector<TripletT> tripletList;
		for(int i = 0; i < n; i++)
		{
			std::map<uint32_t, T> m = t(basis[i]);
			auto iter = basis.begin();
			for(auto& kv: m)
			{
				iter = std::lower_bound(iter, basis.end(), kv.first);
				if(iter == basis.end())
					break;
				auto j = std::distance(basis.begin(), iter);
				{
					tripletList.emplace_back(i, j, kv.second);
				}
			}
		}

		Eigen::SparseMatrix<T> res(n, n);
		res.setFromTriplets(tripletList.begin(), tripletList.end());
		return res;
	}

}

#endif//EDP_CONSTRUCTSPARSEMAT_HPP
