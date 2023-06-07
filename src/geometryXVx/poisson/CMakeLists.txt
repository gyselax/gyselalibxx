# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("poisson_${GEOMETRY_VARIANT}" STATIC
    chargedensitycalculator.cpp
    electricfield.cpp
    nullpoissonsolver.cpp
)

target_compile_features("poisson_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("poisson_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("poisson_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        sll::splines
        vcx::geometry_${GEOMETRY_VARIANT}
        vcx::speciesinfo
)

add_library("vcx::poisson_${GEOMETRY_VARIANT}" ALIAS "poisson_${GEOMETRY_VARIANT}")

endforeach()

target_sources(poisson_xperiod_vx
    PRIVATE
        fftpoissonsolver.cpp
        femperiodicpoissonsolver.cpp
)

target_sources(poisson_xnonperiod_vx
    PRIVATE
        femnonperiodicpoissonsolver.cpp
)
