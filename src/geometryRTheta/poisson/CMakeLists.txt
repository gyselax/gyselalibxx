# SPDX-License-Identifier: MIT

add_library("poisson_RTheta" STATIC
    nullpoissonsolver.cpp
)

target_compile_features("poisson_RTheta"
    PUBLIC
        cxx_std_17
)

target_include_directories("poisson_RTheta"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("poisson_RTheta"
    PUBLIC
        DDC::DDC
        sll::splines
        gslx::geometry_RTheta
        gslx::speciesinfo
)

add_library("gslx::poisson_RTheta" ALIAS "poisson_RTheta")
