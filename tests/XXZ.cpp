#include "XXZ.hpp"

std::map<uint32_t, double> XXZ::operator()(uint32_t col) const
{
    std::map<uint32_t, double> m;
    for(uint32_t i = 0; i < n_; i++)
    {
        uint32_t b1 = (col >> i) & 1U;
        uint32_t b2 = (col >> ((i + 1) % n_)) & 1U;
        int zz = (1 - 2 * static_cast<int>(b1)) * (1 - 2 * static_cast<int>(b2));
        uint64_t x = (1U << i) | (1U << ((i + 1) % (n_)));
        m[col] += J_ * Delta_ * zz;
        m[col ^ x] += J_ * (1 - zz);
    }
    return m;
}
