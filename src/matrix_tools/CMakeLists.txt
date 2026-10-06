
add_library(matrix_tools STATIC
    matrix.cpp
    matrix_dense.cpp
    matrix_banded.cpp
    matrix_corner_block.cpp
    matrix_centre_block.cpp
    matrix_periodic_banded.cpp
    matrix_pds_tridiag.cpp
)

target_compile_features("matrix_tools" PUBLIC cxx_std_20)

target_sources("matrix_tools"
    INTERFACE
        FILE_SET HEADERS
            FILES
                matrix.hpp
                matrix_banded.hpp
                matrix_batch.hpp
                matrix_batch_csr.hpp
                matrix_batch_ell.hpp
                matrix_batch_tridiag.hpp
                matrix_centre_block.hpp
                matrix_corner_block.hpp
                matrix_dense.hpp
                matrix_pds_tridiag.hpp
                matrix_periodic_banded.hpp
                matrix_utils.hpp
)

target_link_libraries(matrix_tools
    PUBLIC
        DDC::core
        DDC::splines
        Kokkos::kokkos
        Ginkgo::ginkgo
        gyselalibxx::utils
        LAPACKE::LAPACKE
)
add_library(gyselalibxx::matrix_tools ALIAS matrix_tools)

install(TARGETS matrix_tools EXPORT gyselalibxxTargets FILE_SET HEADERS)
