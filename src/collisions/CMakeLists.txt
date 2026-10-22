# SPDX-License-Identifier: MIT

add_library("collisions" INTERFACE)

target_compile_features("collisions" INTERFACE cxx_std_20)

target_sources("collisions"
    INTERFACE
        FILE_SET HEADERS
            FILES
                collision_common_configuration.hpp
                collision_operator.hpp
)

target_link_libraries("collisions"
    INTERFACE
        DDC::core
        gyselalibxx::paraconfpp
        gyselalibxx::speciesinfo
        gyselalibxx::utils
        koliop::koliop
        paraconf::paraconf
)

add_library("gyselalibxx::collisions" ALIAS "collisions")

install(TARGETS collisions EXPORT gyselalibxxTargets FILE_SET HEADERS)
