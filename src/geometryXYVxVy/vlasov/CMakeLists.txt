# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("vlasov_xyvxvy_${PRECISION_VARIANT}" STATIC
    mpisplitvlasovsolver.cpp
    splitvlasovsolver.cpp
)

target_sources("vlasov_xyvxvy_${PRECISION_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                ivlasovsolver.hpp
                mpisplitvlasovsolver.hpp
                splitvlasovsolver.hpp
)

target_link_libraries("vlasov_xyvxvy_${PRECISION_VARIANT}"
    PUBLIC
        DDC::core
        gslx::advection
        gslx::geometry_xyvxvy_${PRECISION_VARIANT}
        gslx::mpi_parallelisation
        gslx::speciesinfo
)

add_library("gslx::vlasov_xyvxvy_${PRECISION_VARIANT}" ALIAS "vlasov_xyvxvy_${PRECISION_VARIANT}")

install(TARGETS vlasov_xyvxvy_${PRECISION_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()
