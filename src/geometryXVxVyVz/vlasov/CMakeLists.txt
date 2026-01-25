# SPDX-License-Identifier: MIT

add_library("vlasov_xvxvyvz" STATIC
    mpisplitvlasovsolver.cpp
    mpisplithybridvlasovsolver.cpp
    splitvlasovsolver.cpp
)

target_include_directories("vlasov_xvxvyvz"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("vlasov_xvxvyvz"
    PUBLIC
        DDC::core
        gslx::advection
        gslx::geometry_xvxvyvz
        gslx::mpi_parallelisation
        gslx::speciesinfo
)

add_library("gslx::vlasov_xvxvyvz" ALIAS "vlasov_xvxvyvz")
