cmake_minimum_required(VERSION 3.14)
project(ExactDiagonalization)

################################################################################
# Set options
################################################################################

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_TESTS "Build Tests" ON)

################################################################################
# Set dependencies
################################################################################

find_package(TBB QUIET)
if (NOT TBB_FOUND)
    message(STATUS "TBB is not found. Downloading.")
    include(FetchContent)
    FetchContent_Declare(
        TBB
        GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
        GIT_TAG        v2021.5.0
    )
    FetchContent_Populate(TBB)
    FetchContent_GetProperties(TBB)

    set(TBB_TEST OFF CACHE BOOL " " FORCE)
    add_subdirectory(${tbb_SOURCE_DIR} ${tbb_BINARY_DIR} EXCLUDE_FROM_ALL)
else()
    message(STATUS "TBB found. Use system TBB.")
endif()
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(arpack-ng REQUIRED)

################################################################################
# Set targets
################################################################################

add_library(edlib_dependency INTERFACE)
target_compile_features(edlib_dependency INTERFACE cxx_std_17)
target_compile_options(edlib_dependency INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-Werror>)

add_subdirectory(include)

#export(TARGETS edlib
#       NAMESPACE EDLIB::
#       FILE edlibTargets.cmake)

if (BUILD_EXAMPLES)
	add_executable(tfi_arpack "examples/tfi_arpack.cpp")
	target_link_libraries(tfi_arpack PRIVATE edlib edlib_dependency)
endif()

if (BUILD_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif()
