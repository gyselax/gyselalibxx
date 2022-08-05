# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_VARIANTS_LIST)

add_library("vlasov_${GEOMETRY_VARIANT}" STATIC
    bsl_advection_x.cpp
    bsl_advection_vx.cpp
    nulladvectionvx.cpp
    splitvlasovsolver.cpp
)

target_compile_features("vlasov_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("vlasov_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("vlasov_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        sll::splines
        vcx::geometry_${GEOMETRY_VARIANT}
        vcx::interpolation_${GEOMETRY_VARIANT}
        vcx::speciesinfo
)

add_library("vcx::vlasov_${GEOMETRY_VARIANT}" ALIAS "vlasov_${GEOMETRY_VARIANT}")

endforeach()
