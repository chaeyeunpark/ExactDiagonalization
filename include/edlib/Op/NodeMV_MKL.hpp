#pragma once

#include "exceptions.hpp"

#include <mkl.h>

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

    sparse_matrix_t A_;
    matrix_descr descA_;

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
        sparse_status_t m
            = mkl_sparse_d_create_csr(&A_, SPARSE_INDEX_BASE_ZERO, rows_, dim, ptrB.data(),
                                      ptrB.data() + 1, colIdx.data(), values.data());

        if(m == SPARSE_STATUS_ALLOC_FAILED)
        {
            throw SparseAllocFailed();
        }
        else if(m != SPARSE_STATUS_SUCCESS)
        {
            throw SparseCreateFailed();
        }

        descA_.type = SPARSE_MATRIX_TYPE_GENERAL;

        mkl_sparse_set_mv_hint(A_, SPARSE_OPERATION_NON_TRANSPOSE, descA_, 1024 * 1024);
    }

    [[nodiscard]] size_t row_start() const { return row_start_; }
    [[nodiscard]] std::size_t rows() const { return rows_; }
    [[nodiscard]] std::size_t cols() const { return dim_; }

    void perform_op(const double* x_in, double* y_out) const
    {
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, A_, descA_, x_in, 0.0, y_out);
    }

    NodeMV(const NodeMV& rhs) = default;
    NodeMV(NodeMV&& rhs) noexcept = default;

    NodeMV& operator=(const NodeMV& rhs) = default;
    NodeMV& operator=(NodeMV&& rhs) noexcept = default;

    ~NodeMV() { mkl_sparse_destroy(A_); }
};
} // namespace edlib
