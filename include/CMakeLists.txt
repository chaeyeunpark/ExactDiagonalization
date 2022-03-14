project(ExactDiagonalization-include)

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

add_library(edlib INTERFACE)
target_include_directories(edlib INTERFACE ".")
target_link_libraries(edlib INTERFACE Eigen3::Eigen TBB::tbb ARPACK::ARPACK)

set_target_properties (edlib PROPERTIES EXPORT_NAME edlib)
