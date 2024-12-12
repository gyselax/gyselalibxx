# SPDX-License-Identifier: MIT


add_library(poisson_RTheta INTERFACE)
target_include_directories(poisson_RTheta
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(poisson_RTheta INTERFACE
        DDC::core
        sll::SLL
        gslx::data_types
        gslx::geometry_RTheta
        gslx::speciesinfo
        gslx::utils
        gslx::interpolation_2D_rp
)

add_library("gslx::poisson_RTheta" ALIAS "poisson_RTheta")
