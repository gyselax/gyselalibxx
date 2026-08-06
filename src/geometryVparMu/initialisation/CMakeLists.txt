# SPDX-License-Identifier: MIT

add_library("initialisation_vparmu" STATIC
    maxwellianequilibrium.cpp
    noperturbinitialisation.cpp
)

target_compile_features("initialisation_vparmu" PUBLIC cxx_std_20)

target_sources("initialisation_vparmu"
    INTERFACE
        FILE_SET HEADERS
            FILES
                iequilibrium.hpp
                iinitialisation.hpp
                maxwellianequilibrium.hpp
                noperturbinitialisation.hpp
)

target_link_libraries("initialisation_vparmu"
    PUBLIC
        DDC::core
        gyselalibxx::speciesinfo
        gyselalibxx::geometry_vparmu
        gyselalibxx::utils
)

add_library("gyselalibxx::initialisation_vparmu" ALIAS "initialisation_vparmu")

install(TARGETS initialisation_vparmu EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryVparMu/initialisation)
