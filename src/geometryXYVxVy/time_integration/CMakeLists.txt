# SPDX-License-Identifier: MIT

add_library("time_integration_xyvxvy" STATIC
    predcorr.cpp
)

target_sources("time_integration_xyvxvy"
    INTERFACE
        FILE_SET HEADERS
            FILES
                itimesolver.hpp
                predcorr.hpp
)

target_link_libraries("time_integration_xyvxvy"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::poisson_xy
        gslx::speciesinfo
        gslx::vlasov_xyvxvy
        gslx::utils

)

add_library("gslx::time_integration_xyvxvy" ALIAS "time_integration_xyvxvy")

install(TARGETS time_integration_xyvxvy EXPORT GyselalibxxTargets FILE_SET HEADERS)
