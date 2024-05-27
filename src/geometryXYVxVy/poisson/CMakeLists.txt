# SPDX-License-Identifier: MIT

add_library("poisson_xy" STATIC
    chargedensitycalculator.cpp
    nullqnsolver.cpp
    fftqnsolver.cpp
)

target_compile_features("poisson_xy"
    PUBLIC
        cxx_std_17
)

target_include_directories("poisson_xy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("poisson_xy"
    PUBLIC
        DDC::DDC
        sll::SLL
        gslx::geometry_xyvxvy
        gslx::speciesinfo
        gslx::utils
)

add_library("gslx::poisson_xy" ALIAS "poisson_xy")
