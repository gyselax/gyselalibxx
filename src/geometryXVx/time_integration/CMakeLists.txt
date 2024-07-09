# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("time_integration_${GEOMETRY_VARIANT}" STATIC
    predcorr.cpp
)

target_include_directories("time_integration_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("time_integration_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::boltzmann_${GEOMETRY_VARIANT}
)

add_library("gslx::time_integration_${GEOMETRY_VARIANT}" ALIAS "time_integration_${GEOMETRY_VARIANT}")

endforeach()
