# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("initialization_${GEOMETRY_VARIANT}" STATIC
    maxwellianequilibrium.cpp
)

target_include_directories("initialization_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialization_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        gslx::speciesinfo
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::utils
)

add_library("gslx::initialization_${GEOMETRY_VARIANT}" ALIAS "initialization_${GEOMETRY_VARIANT}")

endforeach()
