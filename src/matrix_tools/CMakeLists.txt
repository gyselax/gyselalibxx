
add_library(matrix_tools STATIC
    matrix.cpp
    matrix_dense.cpp
    matrix_banded.cpp
    matrix_corner_block.cpp
    matrix_centre_block.cpp
    matrix_periodic_banded.cpp
    matrix_pds_tridiag.cpp
)
target_compile_features(matrix_tools PUBLIC cxx_std_17)
target_include_directories(matrix_tools
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
    SYSTEM INTERFACE "${LAPACKE_INCLUDE_DIRS}"
)
target_link_libraries(matrix_tools
    PUBLIC
        DDC::core
        DDC::splines
        Kokkos::kokkos
        gslx::utils
        ${LAPACKE_LIBRARIES}
)
add_library(gslx::matrix_tools ALIAS matrix_tools)
