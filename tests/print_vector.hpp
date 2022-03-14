#pragma once
#include <ostream>
#include <vector>

template<typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    for(const auto& elt : v)
    {
        os << elt << ", ";
    }
    os << "\n";
    return os;
}
