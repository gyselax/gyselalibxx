# SPDX-License-Identifier: MIT

add_library("time_integration_xyvxvy" STATIC
    predcorr.cpp
)

target_compile_features("time_integration_xyvxvy"
    PUBLIC
        cxx_std_17
)

target_include_directories("time_integration_xyvxvy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("time_integration_xyvxvy"
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        vcx::poisson_xy
        vcx::speciesinfo
        vcx::vlasov_xyvxvy
)

add_library("vcx::time_integration_xyvxvy" ALIAS "time_integration_xyvxvy")

