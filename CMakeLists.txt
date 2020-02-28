cmake_minimum_required(VERSION 3.9)

project(ED)
include_directories("${CMAKE_SOURCE_DIR}/Include")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set (CMAKE_CXX_STANDARD 14)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)


set(CMAKE_CXX_FLAGS "-Wall -Wextra -DEIGEN_USE_MKL_ALL")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(EDTFI "Examples/EDTFI.cpp")
target_link_libraries(EDTFI mkl_intel_lp64 mkl_tbb_thread mkl_core tbb pthread m dl Eigen3::Eigen)
