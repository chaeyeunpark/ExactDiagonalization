cmake_minimum_required(VERSION 3.14)

project(ED)
include_directories("${PROJECT_SOURCE_DIR}/externals/spectra/include")
include_directories("${PROJECT_SOURCE_DIR}/externals")

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
        COMMAND            ${PROJECT_SOURCE_DIR}/bin/cpp_files.py --include-examples --include-tests
        WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE    ED_SOURCE_FILES)
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_BINARY};
                             -extra-arg=-std=c++17;
                             --line-filter=${ED_SOURCE_FILES};
    )
endif()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(edlib_dependency INTERFACE)
target_compile_features(edlib_dependency INTERFACE cxx_std_17)
target_compile_options(edlib_dependency INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-Werror>)

add_library(edlib INTERFACE)
target_include_directories(edlib INTERFACE "${PROJECT_SOURCE_DIR}/include")

# MKL is required
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel")
else()
    target_include_directories(edlib SYSTEM INTERFACE "$ENV{MKLROOT}/include")
	target_link_libraries(edlib INTERFACE mkl_intel_lp64 mkl_tbb_thread mkl_core)
endif()

target_link_libraries(edlib INTERFACE tbb pthread m dl)

if (BUILD_EXAMPLES)
	add_executable(tfi "examples/tfi.cpp")
	target_link_libraries(tfi PRIVATE edlib Eigen3::Eigen)
endif()

if (BUILD_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif()
