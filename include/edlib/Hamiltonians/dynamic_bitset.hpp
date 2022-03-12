#pragma once

#include <cstdint>
#include <cstdlib>
#include <vector>

namespace edlib
{
struct dynamic_bitset
{
    std::vector<uint8_t> data_{};

    dynamic_bitset(size_t size, size_t val)
    {
        data_.reserve(size);
        for(size_t i = 0; i < size; i++)
        {
            data_.emplace_back(val & 1U);
            val >>= 1U;
        }
    }

    uint8_t& operator[](size_t i) { return data_[i]; }

    uint8_t operator[](size_t i) const { return data_[i]; }
};
} // namespace edlib
