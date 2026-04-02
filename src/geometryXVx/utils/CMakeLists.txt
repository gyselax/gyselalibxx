# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("utils_${GEOMETRY_VARIANT}" STATIC
    fluid_moments.cpp
)

target_sources("utils_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                fluid_moments.hpp
)

target_link_libraries("utils_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::quadrature
        gslx::utils

)

add_library("gslx::utils_${GEOMETRY_VARIANT}" ALIAS "utils_${GEOMETRY_VARIANT}")

install(TARGETS utils_${GEOMETRY_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()
