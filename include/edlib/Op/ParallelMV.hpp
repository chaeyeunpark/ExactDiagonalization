#pragma once
#include "NodeMV.hpp"

#include <tbb/tbb.h>

#include <cstdlib>

namespace edlib
{
class ParallelMV
{
private:
    uint32_t nDiv_;
    size_t dim_;
    tbb::concurrent_vector<NodeMV> mvs_;

public:
    template<typename ColFunc>
    ParallelMV(const std::size_t dim, ColFunc&& col, int32_t nDiv = -1) : dim_(dim)
    {
        if(nDiv < 0)
        {
            nDiv_ = tbb::this_task_arena::max_concurrency();
        }
        else
        {
            nDiv_ = uint32_t(nDiv);
        }
        tbb::parallel_for(uint32_t(0U), nDiv_, [&](uint32_t idx) {
            mvs_.emplace_back(dim, (dim * idx) / nDiv_, (dim * (idx + 1)) / nDiv_, col);
        });
    }

    void perform_op(const double* x_in, double* y_out) const
    {
        tbb::parallel_for(uint32_t(0U), nDiv_, [&](uint32_t idx) {
            mvs_[idx].perform_op(x_in, y_out + mvs_[idx].row_start());
        });
    }

    [[nodiscard]] size_t rows() const { return dim_; }
    [[nodiscard]] size_t cols() const { return dim_; }
};

} // namespace edlib
