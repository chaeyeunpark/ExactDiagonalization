project(ExactDiagonalization-include)

add_library(edlib INTERFACE)
target_include_directories(edlib INTERFACE ".")
target_link_libraries(edlib INTERFACE Eigen3::Eigen TBB::tbb ARPACK::ARPACK)

set_target_properties (edlib PROPERTIES EXPORT_NAME edlib)
