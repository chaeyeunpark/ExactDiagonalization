cmake_minimum_required(VERSION 3.14)

project(ED)

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

find_package(TBB REQUIRED)
if(NOT TBB_FOUND)
    include(FetchContent)
    message(STATUE "TBB is not found. Downloading.")
    FetchContent_Declare(
        TBB
        GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
        GIT_TAG        v2021.5.0
    )
    set(TBB_BUILD_TESTS OFF)
    FetchContent_MakeAvailable(TBB)
else()
    message(STATUE "TBB found.")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(BLAS REQUIRED)

add_library(edlib_dependency INTERFACE)
target_compile_features(edlib_dependency INTERFACE cxx_std_17)
target_compile_options(edlib_dependency INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-Werror>)

add_library(edlib INTERFACE)
target_include_directories(edlib INTERFACE "${PROJECT_SOURCE_DIR}/include" ${TBB_SOURCE_DIR})
target_link_libraries(edlib INTERFACE "${BLAS_LIBRARIES}" Eigen3::Eigen TBB::tbb)
target_link_options(edlib INTERFACE "${BLAS_LINKER_FLAGS}")


if (BUILD_EXAMPLES)
	add_executable(tfi "examples/tfi.cpp")
	target_link_libraries(tfi PRIVATE edlib)
endif()

if (BUILD_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif()
