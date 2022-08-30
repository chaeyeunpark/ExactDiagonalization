#include "edlib/Basis/BasisJz.hpp"

#include <catch2/catch_all.hpp>

#include <iterator>

TEST_CASE("static check BasisJz iterator", "[BasisJz]")
{
    STATIC_CHECK(std::input_iterator<edlib::BasisJz<uint32_t>::BasisJzIterator>);
}
