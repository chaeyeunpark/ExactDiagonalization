project(ExactDiagonalizaion-Tests)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_library(ed_tests_dependency INTERFACE)
target_link_libraries(ed_tests_dependency INTERFACE Catch2::Catch2 edlib_dependency)

# Required for catch_discover_tests().
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)

# Modify `ctest` to only run the supported subset of tests.
include(CTest)
include(Catch)

add_executable(runner catch_main.cpp
                      test_xxz_gs.cpp
                      test_basisfull.cpp
                      test_basis1d.cpp
                      test_basis2d.cpp
                      test_local_ham.cpp
                      XXZ.cpp
                      utils.cpp)
target_link_libraries(runner PRIVATE edlib Eigen3::Eigen ed_tests_dependency)

catch_discover_tests(runner)

# add_executable(test_perf_basis2d "test_perf_basis2d.cpp")
# target_link_libraries(test_perf_basis2d PRIVATE edlib Eigen3::Eigen)
# 
# add_executable(test_perf_basis1dz2 "test_perf_basis1dz2.cpp")
# target_link_libraries(test_perf_basis1dz2 PRIVATE edlib Eigen3::Eigen)
# 
# add_executable(test_perf_basis2dz2 "test_perf_basis2dz2.cpp")
# target_link_libraries(test_perf_basis2dz2 PRIVATE edlib Eigen3::Eigen)
