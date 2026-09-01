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
        gyselalibxx::poisson_xy_${PRECISION_VARIANT}
        gyselalibxx::speciesinfo
        gyselalibxx::vlasov_xyvxvy_${PRECISION_VARIANT}
        gyselalibxx::utils

)

add_library("gyselalibxx::time_integration_xyvxvy_${PRECISION_VARIANT}" ALIAS "time_integration_xyvxvy_${PRECISION_VARIANT}")

install(TARGETS time_integration_xyvxvy_${PRECISION_VARIANT} EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXYVxVy/time_integration)

endforeach()
