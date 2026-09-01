# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("utils_${GEOMETRY_VARIANT}" STATIC
    fluid_moments.cpp
)

target_compile_features("utils_${GEOMETRY_VARIANT}" PUBLIC cxx_std_20)

target_sources("utils_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                fluid_moments.hpp
)

target_link_libraries("utils_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        gyselalibxx::geometry_${GEOMETRY_VARIANT}
        gyselalibxx::quadrature
        gyselalibxx::utils

)

add_library("gyselalibxx::utils_${GEOMETRY_VARIANT}" ALIAS "utils_${GEOMETRY_VARIANT}")

install(TARGETS utils_${GEOMETRY_VARIANT} EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXVx/utils)

endforeach()
