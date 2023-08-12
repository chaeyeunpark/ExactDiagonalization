#include "edlib/Basis/BasisJz.hpp"

#include <catch2/catch_all.hpp>

#include <iterator>

TEST_CASE("Tests BasisJzIterator") {
    STATIC_CHECK(std::input_iterator<edlib::BasisJz<uint32_t>::BasisJzIterator>);
    STATIC_CHECK(std::forward_iterator<edlib::BasisJz<uint32_t>::BasisJzIterator>);
}
