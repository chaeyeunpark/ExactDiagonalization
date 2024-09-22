cmake_minimum_required(VERSION 3.16)
project(ExactDiagonalization)

################################################################################
# Set options
################################################################################

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

option(EDLIB_BUILD_TESTS "Build Tests" OFF)
option(EDLIB_BUILD_EXAMPLES "Build examples" OFF)
option(ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
option(EDLIB_ENABLE_ARPACK_SOLVER "Enable arpack based solver (ARPACK required)" ON)

################################################################################
# Set dependencies
################################################################################

find_package(TBB QUIET)
if (NOT TBB_FOUND)
    message(STATUS "TBB is not found. Downloading.")
    include(FetchContent)
    FetchContent_Declare(tbb
                         GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
                         GIT_TAG        v2021.5.0
    )
    set(TBB_TEST OFF CACHE BOOL " " FORCE)
    FetchContent_GetProperties(tbb)
    if(NOT tbb_POPULATED)
        FetchContent_Populate(tbb)
        add_subdirectory(${tbb_SOURCE_DIR} ${tbb_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()

    get_target_property(TBB_INC_DIR tbb INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(tbb PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${TBB_INC_DIR}")
else()
    message(STATUS "TBB found. Use system TBB.")
endif()
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

<<<<<<< HEAD
if (EDLIB_ENABLE_ARPACK_SOLVER)
    message(STATUS "ARPACK solver is enabled. Finding ARPACK-NG.")
    find_package(arpack-ng QUIET)
    if (NOT arpack-ng_FOUND)
        find_package(arpackng QUIET)
    endif()

    if ((NOT arpack-ng_FOUND) AND (NOT arpackng_FOUND))
        message(FATAL_ERROR "Arpack-ng is not found. Please provide a path to arpackng by arpackng_ROOT.")
=======
find_package(arpack QUIET)

if (NOT ARPACK_FOUND)
    find_package(ARPACK-NG QUIET)
    find_package(ARPACKNG QUIET)

    if ((NOT ARPACK-NG_FOUND) AND (NOT ARPACKNG_FOUND))
        message(FATAL_ERROR "ARPACK must be available. If ARPACK-NG is installed, please add the path to arpackng-config.cmake in the enviroment variable CMAKE_PREFIX_PATH")
>>>>>>> origin/master
    endif()
endif()

################################################################################
# Set targets
################################################################################

add_library(edlib_dependency INTERFACE)
target_compile_features(edlib_dependency INTERFACE cxx_std_20)
target_compile_options(edlib_dependency INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-Werror>)

# This will include edlib
add_subdirectory(include)
add_library(EDLIB::edlib ALIAS edlib)

if (EDLIB_BUILD_EXAMPLES)
    if(NOT EDLIB_ENABLE_ARPACK_SOLVER)
        message(FATAL_ERROR "Examples can be only built when arpack solver is available.")
    endif()
	add_executable(tfi_arpack "examples/tfi_arpack.cpp")
	target_link_libraries(tfi_arpack PRIVATE edlib edlib_dependency)
endif()

if (EDLIB_BUILD_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif()
