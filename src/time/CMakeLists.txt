# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_VARIANTS_LIST)

add_library("time_${GEOMETRY_VARIANT}" STATIC
    predcorr.cpp
    rk2_solver.cpp
)

target_compile_features("time_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("time_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("time_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        vcx::poisson_${GEOMETRY_VARIANT}
        vcx::speciesinfo
        vcx::boltzmann_${GEOMETRY_VARIANT}
)

add_library("vcx::time_${GEOMETRY_VARIANT}" ALIAS "time_${GEOMETRY_VARIANT}")

endforeach()
