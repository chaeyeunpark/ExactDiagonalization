project(ExactDiagonalizaion-Tests)

enable_testing() 

add_executable(test_xxz_gs "test_xxz_gs.cpp")
target_link_libraries(test_xxz_gs PRIVATE edlib Eigen3::Eigen)
add_test(    
	NAME test_xxz_gs
	COMMAND test_xxz_gs
)

add_executable(test_basisfull "test_basisfull.cpp")
target_link_libraries(test_basisfull PRIVATE edlib Eigen3::Eigen)
add_test(    
	NAME test_basisfull
	COMMAND test_basisfull
)

add_executable(test_basis1d "test_basis1d.cpp")
target_link_libraries(test_basis1d PRIVATE edlib Eigen3::Eigen)
add_test(    
	NAME test_basis1d
	COMMAND test_basis1d
)

add_executable(test_basis2d "test_basis2d.cpp")
target_link_libraries(test_basis2d PRIVATE edlib Eigen3::Eigen)
add_test(    
	NAME test_basis2d
	COMMAND test_basis2d
)

add_executable(test_local_ham "test_local_ham.cpp")
target_link_libraries(test_local_ham PRIVATE edlib Eigen3::Eigen)
add_test(    
	NAME test_local_ham
	COMMAND test_local_ham
)

add_executable(test_perf_basis2d "test_perf_basis2d.cpp")
target_link_libraries(test_perf_basis2d PRIVATE edlib Eigen3::Eigen)

add_executable(test_perf_basis1dz2 "test_perf_basis1dz2.cpp")
target_link_libraries(test_perf_basis1dz2 PRIVATE edlib Eigen3::Eigen)

add_executable(test_perf_basis2dz2 "test_perf_basis2dz2.cpp")
target_link_libraries(test_perf_basis2dz2 PRIVATE edlib Eigen3::Eigen)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
	DEPENDS test_xxz_gs test_basisfull test_basis1d test_basis2d test_local_ham)


