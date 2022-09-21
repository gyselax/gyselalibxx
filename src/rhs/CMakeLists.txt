# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_VARIANTS_LIST)

add_library("rhs_${GEOMETRY_VARIANT}" STATIC
    kinetic_source.cpp
)

target_compile_features("rhs_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("rhs_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("rhs_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        vcx::speciesinfo
        vcx::geometry_${GEOMETRY_VARIANT}
)

add_library("vcx::rhs_${GEOMETRY_VARIANT}" ALIAS "rhs_${GEOMETRY_VARIANT}")

endforeach()
