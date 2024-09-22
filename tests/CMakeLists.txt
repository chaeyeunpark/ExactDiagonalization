project(ExactDiagonalization-Tests)

Include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.1.0
)

FetchContent_MakeAvailable(Catch2)
get_target_property(Catch2_INC_DIR Catch2 INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(Catch2 PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Catch2_INC_DIR}")
# Required for catch_discover_tests().
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

FetchContent_Declare(
    Spectra
    GIT_REPOSITORY https://github.com/yixuan/spectra.git
    GIT_TAG v1.0.0
)
set(BUILD_TESTS OFF CACHE INTERNAL "Turn off tests for spectra")
FetchContent_MakeAvailable(Spectra)
get_target_property(Spectra_INC_DIR Spectra INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(Spectra PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Spectra_INC_DIR}")

if(ENABLE_CLANG_TIDY)
    if(NOT DEFINED CLANG_TIDY_BINARY)
        set(CLANG_TIDY_BINARY clang-tidy-13)
    endif()
    message(STATUS "Using CLANG_TIDY_BINARY=${CLANG_TIDY_BINARY}")
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_BINARY};
                             -extra-arg=-std=c++20;
    )
endif()

add_library(ed_tests_dependency INTERFACE)
target_link_libraries(ed_tests_dependency INTERFACE Catch2::Catch2 Spectra edlib_dependency)

# Modify `ctest` to only run the supported subset of tests.
include(CTest)
include(Catch)

add_executable(runner test_arpack_solver.cpp
                      test_basisfull.cpp
                      test_basis1d.cpp
                      test_basis2d.cpp
                      test_basisjz.cpp
                      test_local_ham.cpp
                      test_op.cpp
                      test_transform_basis.cpp
                      test_xxz_gs.cpp
                      utils.cpp
                      XXZ.cpp)

target_link_libraries(runner PRIVATE edlib Catch2::Catch2WithMain Eigen3::Eigen ed_tests_dependency)

if(EDLIB_ENABLE_ARPACK_SOLVER)
    target_sources(runner PRIVATE test_arpack_solver.cpp)
    target_link_libraries(runner PRIVATE ARPACK::ARPACK)
endif()


catch_discover_tests(runner)

# add_executable(test_perf_basis2d "test_perf_basis2d.cpp")
# target_link_libraries(test_perf_basis2d PRIVATE edlib Eigen3::Eigen)
# 
# add_executable(test_perf_basis1dz2 "test_perf_basis1dz2.cpp")
# target_link_libraries(test_perf_basis1dz2 PRIVATE edlib Eigen3::Eigen)
# 
# add_executable(test_perf_basis2dz2 "test_perf_basis2dz2.cpp")
# target_link_libraries(test_perf_basis2dz2 PRIVATE edlib Eigen3::Eigen)
