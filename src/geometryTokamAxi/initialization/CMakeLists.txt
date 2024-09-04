# SPDX-License-Identifier: MIT

add_library("initialization_tokamaxi" STATIC
    maxwellianequilibrium.cpp
)

target_include_directories("initialization_tokamaxi"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialization_tokamaxi"
    PUBLIC
        DDC::DDC
        gslx::speciesinfo
        gslx::geometry_tokamaxi
)

add_library("gslx::initialization_tokamaxi" ALIAS "initialization_tokamaxi")
