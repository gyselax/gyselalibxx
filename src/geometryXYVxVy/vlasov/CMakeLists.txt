# SPDX-License-Identifier: MIT

add_library("vlasov_xyvxvy" STATIC
    mpisplitvlasovsolver.cpp
    splitvlasovsolver.cpp
)

target_sources("vlasov_xyvxvy"
    INTERFACE
        FILE_SET HEADERS
            FILES
                ivlasovsolver.hpp
                mpisplitvlasovsolver.hpp
                splitvlasovsolver.hpp
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

install(TARGETS vlasov_xyvxvy EXPORT GyselalibxxTargets FILE_SET HEADERS)
