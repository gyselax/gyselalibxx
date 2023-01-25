# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_VARIANTS_LIST)
    
add_library("utils_${GEOMETRY_VARIANT}" STATIC
    fluid_moments.cpp
)

target_compile_features("utils_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("utils_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("utils_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        vcx::geometry_${GEOMETRY_VARIANT}
        vcx::quadrature
)

add_library("vcx::utils_${GEOMETRY_VARIANT}" ALIAS "utils_${GEOMETRY_VARIANT}")

endforeach()
