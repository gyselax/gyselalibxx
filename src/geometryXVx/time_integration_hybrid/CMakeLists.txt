# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("time_integration_hybrid_${GEOMETRY_VARIANT}" STATIC
    predcorr_hybrid.cpp
)

target_include_directories("time_integration_hybrid_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("time_integration_hybrid_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::initialization_${GEOMETRY_VARIANT}
        gslx::fluidsolver_${GEOMETRY_VARIANT}
        gslx::fluidrhs_${GEOMETRY_VARIANT}
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::speciesinfo
)

add_library("gslx::time_integration_hybrid_${GEOMETRY_VARIANT}" ALIAS "time_integration_hybrid_${GEOMETRY_VARIANT}")

endforeach()
