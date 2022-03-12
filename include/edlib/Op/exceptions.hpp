#pragma once
#include <exception>

namespace edlib
{
class SparseAllocFailed : public std::exception
{
    const char* what() const noexcept { return "Allocation failed for a sparse matrix."; }
};
class SparseCreateFailed : public std::exception
{
    const char* what() const noexcept { return "Failed to create a sparse matrix"; }
};
}
