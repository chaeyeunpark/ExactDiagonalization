#pragma once
#include <exception>

namespace edlib
{
class SparseAllocFailed : public std::exception
{
    [[nodiscard]] const char* what() const noexcept override { return "Allocation failed for a sparse matrix."; }
};
class SparseCreateFailed : public std::exception
{
    [[nodiscard]] const char* what() const noexcept override { return "Failed to create a sparse matrix"; }
};
} // namespace edlib
