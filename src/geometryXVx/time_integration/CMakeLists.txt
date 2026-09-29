# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("time_integration_${GEOMETRY_VARIANT}" STATIC
    predcorr.cpp
)

target_compile_features("time_integration_${GEOMETRY_VARIANT}" PUBLIC cxx_std_20)

target_sources("time_integration_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                itimesolver.hpp
                predcorr.hpp
)

target_link_libraries("time_integration_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gyselalibxx::poisson_${GEOMETRY_VARIANT}
        gyselalibxx::speciesinfo
        gyselalibxx::boltzmann_${GEOMETRY_VARIANT}
        gyselalibxx::utils

)

add_library("gyselalibxx::time_integration_${GEOMETRY_VARIANT}" ALIAS "time_integration_${GEOMETRY_VARIANT}")

install(TARGETS time_integration_${GEOMETRY_VARIANT} EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXVx/time_integration)

endforeach()
