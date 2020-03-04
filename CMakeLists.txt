cmake_minimum_required(VERSION 3.9)

project(ED)
include_directories("${CMAKE_SOURCE_DIR}/Include")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set (CMAKE_CXX_STANDARD 14)
enable_testing() 

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -DEIGEN_USE_MKL_ALL")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

add_executable(EDTFI "Examples/EDTFI.cpp")
target_link_libraries(EDTFI mkl_intel_lp64 mkl_tbb_thread mkl_core tbb pthread m dl Eigen3::Eigen)

add_executable(TestBasis "Tests/TestBasis.cpp")
target_link_libraries(TestBasis mkl_intel_lp64 mkl_tbb_thread mkl_core tbb pthread m dl Eigen3::Eigen)
add_test(    
    NAME TestBasis
    COMMAND TestBasis
)

add_executable(TestBasisGS "Tests/TestBasisGS.cpp")
target_link_libraries(TestBasisGS mkl_intel_lp64 mkl_tbb_thread mkl_core tbb pthread m dl Eigen3::Eigen)
add_test(    
    NAME TestBasisGS
    COMMAND TestBasisGS
)

add_executable(TestLocalHam "Tests/TestLocalHam.cpp")
target_link_libraries(TestLocalHam mkl_intel_lp64 mkl_tbb_thread mkl_core tbb pthread m dl Eigen3::Eigen)
add_test(    
    NAME TestLocalHam
    COMMAND TestLocalHam
)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}    
    DEPENDS TestBasis TestBasisGS TestLocalHam)
