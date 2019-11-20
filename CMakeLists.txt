cmake_minimum_required(VERSION 3.9)

project(ED)
include_directories("${CMAKE_SOURCE_DIR}/Include")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-int-in-bool-context -Wno-deprecated-copy -DEIGEN_USE_MKL_ALL")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

link_directories("~/local/lib")
include_directories(PUBLIC ${CMAKE_INCLUDE_PATH})
include_directories(PUBLIC ${EIGEN3_PATH})


add_executable(TestLocalHam "Tests/TestLocalHam.cpp")
target_link_libraries(TestLocalHam mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread dl)

add_executable(TestBasisGS "Tests/TestBasisGS.cpp")
target_link_libraries(TestBasisGS mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread dl)
