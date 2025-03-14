# SPDX-License-Identifier: MIT

add_library("vlasov_xyvxvy" STATIC
    mpisplitvlasovsolver.cpp
    splitvlasovsolver.cpp
)

target_include_directories("vlasov_xyvxvy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("vlasov_xyvxvy"
    PUBLIC
        DDC::core
        gslx::advection
        gslx::geometry_xyvxvy
        gslx::mpi_parallelisation
        gslx::speciesinfo
)

add_library("gslx::vlasov_xyvxvy" ALIAS "vlasov_xyvxvy")
