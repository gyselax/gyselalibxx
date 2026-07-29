# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("time_integration_xyvxvy_${PRECISION_VARIANT}" STATIC
    predcorr.cpp
)

target_compile_features("time_integration_xyvxvy_${PRECISION_VARIANT}" PUBLIC cxx_std_20)

target_sources("time_integration_xyvxvy_${PRECISION_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                itimesolver.hpp
                predcorr.hpp
)

target_link_libraries("time_integration_xyvxvy_${PRECISION_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::poisson_xy_${PRECISION_VARIANT}
        gslx::speciesinfo
        gslx::vlasov_xyvxvy_${PRECISION_VARIANT}
        gslx::utils

)

add_library("gslx::time_integration_xyvxvy_${PRECISION_VARIANT}" ALIAS "time_integration_xyvxvy_${PRECISION_VARIANT}")

install(TARGETS time_integration_xyvxvy_${PRECISION_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXYVxVy/time_integration)

endforeach()
