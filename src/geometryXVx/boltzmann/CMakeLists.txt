# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("boltzmann_${GEOMETRY_VARIANT}" STATIC
    splitvlasovsolver.cpp
    splitrighthandsidesolver.cpp
)

target_compile_features("boltzmann_${GEOMETRY_VARIANT}" PUBLIC cxx_std_20)

target_sources("boltzmann_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                iboltzmannsolver.hpp
                splitrighthandsidesolver.hpp
                splitvlasovsolver.hpp
)

target_link_libraries("boltzmann_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        gyselalibxx::interpolation
        gyselalibxx::speciesinfo
        gyselalibxx::geometry_${GEOMETRY_VARIANT}
        gyselalibxx::rhs_${GEOMETRY_VARIANT}
        gyselalibxx::advection
)

add_library("gyselalibxx::boltzmann_${GEOMETRY_VARIANT}" ALIAS "boltzmann_${GEOMETRY_VARIANT}")

install(TARGETS boltzmann_${GEOMETRY_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXVx/boltzmann)

endforeach()
