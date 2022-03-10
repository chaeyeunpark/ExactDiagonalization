cmake_minimum_required(VERSION 3.9)

project(ED)
include_directories("${PROJECT_SOURCE_DIR}/externals/spectra/include")
include_directories("${PROJECT_SOURCE_DIR}/externals")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CATCH_URL "https://github.com/catchorg/Catch2/releases/download/v2.11.1/catch.hpp")
set(CATCH_DOWNLOAD_PATH "${PROJECT_SOURCE_DIR}/externals/catch.hpp")
if (NOT EXISTS "${CATCH_DOWNLOAD_PATH}")
	file(DOWNLOAD "${CATCH_URL}" "${CATCH_DOWNLOAD_PATH}")
endif()

option(BUILD_TESTS "Build Tests" ON)

set (CMAKE_CXX_STANDARD 17)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_library(edlib INTERFACE)
target_include_directories(edlib INTERFACE "${PROJECT_SOURCE_DIR}/include")

# MKL is required
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel")
else()
	include_directories("$ENV{MKLROOT}/include")
	target_link_libraries(edlib INTERFACE mkl_intel_lp64 mkl_tbb_thread mkl_core)
endif()

target_link_libraries(edlib INTERFACE tbb pthread m dl)

if (BUILD_EXAMPLES)
	add_executable(tfi "examples/tfi.cpp")
	target_link_libraries(tfi PRIVATE edlib Eigen3::Eigen)
endif()

if (BUILD_TESTS)
	add_subdirectory(tests)
endif()
