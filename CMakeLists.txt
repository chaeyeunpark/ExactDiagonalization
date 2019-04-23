cmake_minimum_required(VERSION 3.9)

project(ED)
include_directories(${CMAKE_SOURCE_DIR})

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-int-in-bool-context ")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#set (CMAKE_CXX_FLAGS_RELEASE "-fopenmp")
link_directories("$HOME/local/lib")
include_directories(PUBLIC "~/local/include/")

add_executable(TestBasis "Tests/TestBasis.cpp")
#target_link_libraries(RBMExact )

add_executable(EDTFI "EDTFI.cpp")
target_link_libraries(EDTFI mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread dl)

add_executable(TIXXZED "TIXXZED.cpp")
target_link_libraries(TIXXZED mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread dl)

add_executable(TIXYZED "TIXYZED.cpp")
target_link_libraries(TIXYZED mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread dl)

add_executable(ED_XXXJ1J2 "ED_XXXJ1J2.cpp")
target_link_libraries(ED_XXXJ1J2 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread dl)
