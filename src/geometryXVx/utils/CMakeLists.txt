# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)
    
add_library("utils_${GEOMETRY_VARIANT}" STATIC
    fluid_moments.cpp
)

target_include_directories("utils_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("utils_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::quadrature
)

add_library("gslx::utils_${GEOMETRY_VARIANT}" ALIAS "utils_${GEOMETRY_VARIANT}")

endforeach()
