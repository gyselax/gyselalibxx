# SPDX-License-Identifier: MIT

add_library("pde_solvers" INTERFACE)

target_compile_features("pde_solvers" INTERFACE cxx_std_20)

target_sources("pde_solvers"
    INTERFACE
        FILE_SET HEADERS
            FILES
                fem_1d_poisson_solver.hpp
                fft_poisson_solver.hpp
                ipoisson_solver.hpp
                polarpoissonlikeassembler.hpp
                polarpoissonlikesolver.hpp
)

target_link_libraries("pde_solvers"
    INTERFACE
        DDC::core
        DDC::fft
        gslx::coord_transformations
        gslx::data_types
        gslx::matrix_tools
        gslx::utils
)

add_library("gslx::pde_solvers" ALIAS "pde_solvers")

install(TARGETS pde_solvers EXPORT GyselalibxxTargets FILE_SET HEADERS)
