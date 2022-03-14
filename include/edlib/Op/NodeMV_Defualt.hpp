#pragma once

#include "exceptions.hpp"

#include <algorithm>
#include <cstdlib>
#include <utility>
#include <vector>

namespace edlib
{
class NodeMV
{
private:
    size_t dim_;
    size_t row_start_;
    size_t rows_;

    std::vector<double> values;
    std::vector<int> colIdx;
    std::vector<int> ptrB;

public:
    using Scalar = double;

    template<class ColFunc>
    NodeMV(size_t dim, size_t row_start, size_t row_end, ColFunc&& col)
        : dim_{dim}, row_start_{row_start}, rows_{row_end - row_start}
    {
        auto get_first = [](const std::pair<const std::size_t, double>& p) {
            return p.first;
        };
        auto get_second = [](const std::pair<const std::size_t, double>& p) {
            return p.second;
        };

        ptrB.resize(rows_ + 1);
        ptrB[0] = 0;
        for(std::size_t i = 0; i < rows_; i++)
        {
            auto rr = col.getCol(i + row_start);
            ptrB[i + 1] = ptrB[i] + rr.size();
            std::transform(rr.begin(), rr.end(), back_inserter(colIdx), get_first);
            std::transform(rr.begin(), rr.end(), back_inserter(values), get_second);
        }
    }

    [[nodiscard]] size_t row_start() const { return row_start_; }
    [[nodiscard]] size_t rows() const { return rows_; }
    [[nodiscard]] size_t cols() const { return dim_; }

    void perform_op(const double* x_in, double* y_out) const
    {
        std::fill_n(y_out, rows_, 0.0);
        for(size_t j = 0; j < rows_; ++j)
        {
            for(int p = ptrB[j]; p < ptrB[j + 1]; ++p)
            {
                y_out[j] += values[p] * x_in[colIdx[p]];
            }
        }
    }
};
} // namespace edlib
