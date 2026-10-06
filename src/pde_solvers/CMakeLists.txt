# SPDX-License-Identifier: MIT

add_library("pde_solvers" INTERFACE)

target_compile_features("pde_solvers" INTERFACE cxx_std_20)

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
        gyselalibxx::coord_transformations
        gyselalibxx::data_types
        gyselalibxx::matrix_tools
        gyselalibxx::utils
)

add_library("gyselalibxx::pde_solvers" ALIAS "pde_solvers")

install(TARGETS pde_solvers EXPORT gyselalibxxTargets FILE_SET HEADERS)
