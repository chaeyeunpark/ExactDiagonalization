#ifndef ED_PARALLELMV_HPP
#define ED_PARALLELMV_HPP
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include "NodeMV.hpp"
class ParallelMV
{
private:
	uint32_t nDiv_;
	std::size_t dim_;
	std::vector<std::unique_ptr<NodeMV>> mvs_;

public:
	template<typename ColFunc>
	ParallelMV(const std::size_t dim, ColFunc&& col, int32_t nDiv = -1)
		: dim_(dim)
	{
		if(nDiv < 0)
			nDiv_ = tbb::task_scheduler_init::default_num_threads();
		else
			nDiv_ = uint32_t(nDiv);
		mvs_.resize(nDiv_);
		tbb::parallel_for(uint32_t(0u), nDiv_, 
			[&](uint32_t idx)
		{
			mvs_[idx] = std::make_unique<NodeMV>(dim, (dim*idx)/nDiv_, (dim*(idx+1))/nDiv_, col);
		});
	}

	void perform_op(const double* x_in, double* y_out) const
	{
		tbb::parallel_for(uint32_t(0u), nDiv_, 
			[&](uint32_t idx)
		{
			mvs_[idx]->perform_op(x_in, y_out + (dim_*idx)/nDiv_);
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
