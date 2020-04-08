cmake_minimum_required(VERSION 3.9)

project(ED)
include_directories("${CMAKE_SOURCE_DIR}/Include" ${Spectra_Root}/include)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set (CMAKE_CXX_STANDARD 14)
enable_testing() 

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -DEIGEN_USE_MKL_ALL")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	include_directories("$ENV{MKLROOT}/include")
	link_libraries(mkl_intel_lp64 mkl_tbb_thread mkl_core)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	include_directories("$ENV{MKLROOT}/include")
	link_libraries(mkl_intel_lp64 mkl_tbb_thread mkl_core)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel")
	link_libraries(pthread m dl)
endif()

link_libraries(tbb pthread m dl)

add_executable(EDTFI "Examples/EDTFI.cpp")
target_link_libraries(EDTFI Eigen3::Eigen)

add_executable(TestBasis2D "Tests/TestBasis2D.cpp")
target_link_libraries(TestBasis2D Eigen3::Eigen)

add_executable(TestBasis "Tests/TestBasis.cpp")
target_link_libraries(TestBasis Eigen3::Eigen)
add_test(    
    NAME TestBasis
    COMMAND TestBasis
)

add_executable(TestBasisGS "Tests/TestBasisGS.cpp")
target_link_libraries(TestBasisGS Eigen3::Eigen)
add_test(    
    NAME TestBasisGS
    COMMAND TestBasisGS
)

add_executable(TestLocalHam "Tests/TestLocalHam.cpp")
target_link_libraries(TestLocalHam Eigen3::Eigen)
add_test(    
    NAME TestLocalHam
    COMMAND TestLocalHam
)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}    
    DEPENDS TestBasis TestBasisGS TestLocalHam)
