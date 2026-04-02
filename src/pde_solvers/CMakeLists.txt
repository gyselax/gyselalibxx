# SPDX-License-Identifier: MIT

add_library("pde_solvers" INTERFACE)

target_sources("pde_solvers"
    INTERFACE
        FILE_SET HEADERS
            FILES
                fem_1d_poisson_solver.hpp
                fft_poisson_solver.hpp
                gmg_polar_poisson_like_solver.hpp
                ipoisson_solver.hpp
                ipolar_poisson_like_solver.hpp
                polar_spline_fem_poisson_like_assembler.hpp
                polar_spline_fem_poisson_like_solver.hpp
)

target_link_libraries("pde_solvers"
    INTERFACE
        DDC::core
        DDC::fft
        GMGPolar::GMGPolarLib
        gslx::coord_transformations
        gslx::data_types
        gslx::matrix_tools
        gslx::utils
)

add_library("gslx::pde_solvers" ALIAS "pde_solvers")

install(TARGETS pde_solvers EXPORT GyselalibxxTargets FILE_SET HEADERS)
