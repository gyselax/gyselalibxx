# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("vlasov_xyvxvy_${PRECISION_VARIANT}" STATIC
    mpisplitvlasovsolver.cpp
    splitvlasovsolver.cpp
)

target_include_directories("vlasov_xyvxvy_${PRECISION_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
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

endforeach()
