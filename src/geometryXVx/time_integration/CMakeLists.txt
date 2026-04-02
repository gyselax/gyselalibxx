# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("time_integration_${GEOMETRY_VARIANT}" STATIC
    predcorr.cpp
)

target_sources("time_integration_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                itimesolver.hpp
                predcorr.hpp
)

target_link_libraries("time_integration_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::utils

)

add_library("gslx::time_integration_${GEOMETRY_VARIANT}" ALIAS "time_integration_${GEOMETRY_VARIANT}")

install(TARGETS time_integration_${GEOMETRY_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()
