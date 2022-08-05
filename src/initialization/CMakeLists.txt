# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_VARIANTS_LIST)

add_library("initialization_${GEOMETRY_VARIANT}" STATIC
    maxwellianequilibrium.cpp
    bumpontailequilibrium.cpp
    singlemodeperturbinitialization.cpp
)

target_compile_features("initialization_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("initialization_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialization_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        vcx::speciesinfo
        vcx::geometry_${GEOMETRY_VARIANT}
)

add_library("vcx::initialization_${GEOMETRY_VARIANT}" ALIAS "initialization_${GEOMETRY_VARIANT}")

endforeach()
