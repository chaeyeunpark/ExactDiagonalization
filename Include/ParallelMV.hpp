#ifndef ED_PARALLELMV_HPP
#define ED_PARALLELMV_HPP
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include "NodeMV.hpp"
class ParallelMV
{
private:
	unsigned int nThreads_;
	std::size_t dim_;
	std::vector<std::unique_ptr<NodeMV>> mvs_;
public:
	template<typename ColFunc>
	ParallelMV(const std::size_t dim, ColFunc&& col)
		: dim_(dim)
	{
		nThreads_ = tbb::task_scheduler_init::default_num_threads();
		mvs_.resize(nThreads_);
		tbb::parallel_for(0u, nThreads_, [&](unsigned int idx)
		{
			mvs_[idx] = std::make_unique<NodeMV>(dim, (dim*idx)/nThreads_, (dim*(idx+1))/nThreads_, col);
		});
	}

	void perform_op(const double* x_in, double* y_out) const
	{
		tbb::parallel_for(0u, nThreads_, [&](unsigned int idx)
		{
			mvs_[idx]->perform_op(x_in + (dim_*idx)/nThreads_, y_out + (dim_*idx)/nThreads_);
		});	
	}

	std::size_t rows() const
	{
		return dim_;
	}
	std::size_t cols() const
	{
		return dim_;
	}

};
#endif//ED_PARALLELMV_HPP
