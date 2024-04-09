# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("fluidsolver_${GEOMETRY_VARIANT}" STATIC
    nullfluidsolver.cpp
)

target_compile_features("fluidsolver_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("fluidsolver_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("fluidsolver_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        sll::SLL
        gslx::advection
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::interpolation
        gslx::quadrature
        gslx::speciesinfo
)

add_library("gslx::fluidsolver_${GEOMETRY_VARIANT}" ALIAS "fluidsolver_${GEOMETRY_VARIANT}")

endforeach()
