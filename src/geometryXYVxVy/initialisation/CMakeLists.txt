# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("initialisation_xyvxvy_${PRECISION_VARIANT}" STATIC
    maxwellianequilibrium.cpp
    singlemodeperturbinitialisation.cpp
)

target_compile_features("initialisation_xyvxvy_${PRECISION_VARIANT}" PUBLIC cxx_std_20)

target_sources("initialisation_xyvxvy_${PRECISION_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                iequilibrium.hpp
                iinitialisation.hpp
                maxwellianequilibrium.hpp
                singlemodeperturbinitialisation.hpp
)

target_link_libraries("initialisation_xyvxvy_${PRECISION_VARIANT}"
    PUBLIC
        DDC::core
        gyselalibxx::speciesinfo
        gyselalibxx::geometry_xyvxvy_${PRECISION_VARIANT}
        gyselalibxx::utils
)

add_library("gyselalibxx::initialisation_xyvxvy_${PRECISION_VARIANT}" ALIAS "initialisation_xyvxvy_${PRECISION_VARIANT}")

install(TARGETS initialisation_xyvxvy_${PRECISION_VARIANT} EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXYVxVy/initialisation)

endforeach()
