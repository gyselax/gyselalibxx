# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("fluidinitialization_${GEOMETRY_VARIANT}" STATIC
    constantfluidinitialization.cpp
)

target_include_directories("fluidinitialization_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("fluidinitialization_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        gslx::geometry_${GEOMETRY_VARIANT}
)

add_library("gslx::fluidinitialization_${GEOMETRY_VARIANT}" ALIAS "fluidinitialization_${GEOMETRY_VARIANT}")

endforeach()
