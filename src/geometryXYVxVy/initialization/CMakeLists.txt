# SPDX-License-Identifier: MIT

add_library("initialization_xyvxvy" STATIC
    maxwellianequilibrium.cpp
    singlemodeperturbinitialization.cpp
)

target_compile_features("initialization_xyvxvy"
    PUBLIC
        cxx_std_17
)

target_include_directories("initialization_xyvxvy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialization_xyvxvy"
    PUBLIC
        DDC::DDC
        vcx::speciesinfo
        vcx::geometry_xyvxvy
        vcx::utils
)

add_library("vcx::initialization_xyvxvy" ALIAS "initialization_xyvxvy")
