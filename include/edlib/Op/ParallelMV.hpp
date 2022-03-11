#pragma once
#include "NodeMV.hpp"
#include <tbb/tbb.h>

namespace edlib
{
class ParallelMV
{
private:
    uint32_t nDiv_;
    std::size_t dim_;
    std::vector<std::unique_ptr<NodeMV>> mvs_;

public:
    template<typename ColFunc>
    ParallelMV(const std::size_t dim, ColFunc&& col, int32_t nDiv = -1) : dim_(dim)
    {
        if(nDiv < 0)
            nDiv_ = tbb::this_task_arena::max_concurrency();
        else
            nDiv_ = uint32_t(nDiv);
        mvs_.resize(nDiv_);
        tbb::parallel_for(uint32_t(0u), nDiv_, [&](uint32_t idx) {
            mvs_[idx] = std::make_unique<NodeMV>(dim, (dim * idx) / nDiv_,
                                                 (dim * (idx + 1)) / nDiv_, col);
        });
    }

    void perform_op(const double* x_in, double* y_out) const
    {
        tbb::parallel_for(uint32_t(0u), nDiv_, [&](uint32_t idx) {
            mvs_[idx]->perform_op(x_in, y_out + (dim_ * idx) / nDiv_);
        });
    }

    std::size_t rows() const { return dim_; }
    std::size_t cols() const { return dim_; }
};

} // namespace edlib
