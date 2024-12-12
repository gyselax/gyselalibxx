# SPDX-License-Identifier: MIT

add_library("collisions" STATIC
    koliop_interface.cpp
)

target_include_directories("collisions"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("collisions"
    PUBLIC
        DDC::core
        gslx::paraconfpp
        gslx::speciesinfo
        gslx::utils
        koliop::koliop
        paraconf::paraconf
)

add_library("gslx::collisions" ALIAS "collisions")
