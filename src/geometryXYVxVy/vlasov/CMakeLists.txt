# SPDX-License-Identifier: MIT

add_library("vlasov_xyvxvy" STATIC
    splitvlasovsolver.cpp
)

target_compile_features("vlasov_xyvxvy"
    PUBLIC
        cxx_std_17
)

target_include_directories("vlasov_xyvxvy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("vlasov_xyvxvy"
    PUBLIC
        DDC::DDC
        sll::SLL
        gslx::speciesinfo
        gslx::geometry_xyvxvy
        gslx::advection
)

add_library("gslx::vlasov_xyvxvy" ALIAS "vlasov_xyvxvy")
