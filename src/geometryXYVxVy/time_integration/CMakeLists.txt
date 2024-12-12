# SPDX-License-Identifier: MIT

add_library("time_integration_xyvxvy" STATIC
    predcorr.cpp
)

target_include_directories("time_integration_xyvxvy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("time_integration_xyvxvy"
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        gslx::poisson_xy
        gslx::speciesinfo
        gslx::vlasov_xyvxvy
        gslx::utils

)

add_library("gslx::time_integration_xyvxvy" ALIAS "time_integration_xyvxvy")

