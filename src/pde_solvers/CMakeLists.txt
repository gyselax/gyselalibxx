# SPDX-License-Identifier: MIT

add_library("pde_solvers" INTERFACE)

target_include_directories("pde_solvers"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
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
