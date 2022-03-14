cmake_minimum_required(VERSION 3.14)
project(ExactDiagonalization)

################################################################################
# Set options
################################################################################

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_TESTS "Build Tests" ON)

if(ENABLE_CLANG_TIDY)
    if(NOT DEFINED CLANG_TIDY_BINARY)
        set(CLANG_TIDY_BINARY clang-tidy-13)
    endif()
    message(STATUS "Using CLANG_TIDY_BINARY=${CLANG_TIDY_BINARY}")
    execute_process(
        COMMAND            ./bin/cpp_files.py --include-examples --include-tests
        WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE    ED_SOURCE_FILES)
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_BINARY};
                             -extra-arg=-std=c++17;
                             --line-filter=${ED_SOURCE_FILES};
    )
endif()

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
    get_target_property(TBB_IID INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(TBB PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${TBB_IID}")
    add_subdirectory(${tbb_SOURCE_DIR} ${tbb_BINARY_DIR})
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

add_library(edlib INTERFACE)
target_include_directories(edlib INTERFACE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(edlib INTERFACE Eigen3::Eigen TBB::tbb ARPACK::ARPACK)

set_target_properties (edlib PROPERTIES EXPORT_NAME edlib)

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
