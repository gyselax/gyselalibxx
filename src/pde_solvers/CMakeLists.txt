# SPDX-License-Identifier: MIT

add_library("pde_solvers" INTERFACE)

target_compile_features("pde_solvers"
    INTERFACE
        cxx_std_17
)

target_include_directories("pde_solvers"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("pde_solvers"
    INTERFACE
        DDC::DDC
        sll::SLL
        gslx::data_types
        gslx::utils
)

add_library("gslx::pde_solvers" ALIAS "pde_solvers")
