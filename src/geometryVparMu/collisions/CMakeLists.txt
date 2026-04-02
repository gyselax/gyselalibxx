# SPDX-License-Identifier: MIT

add_library("geometry_collisions_vparmu" INTERFACE)

target_sources("geometry_collisions_vparmu"
    INTERFACE
        FILE_SET HEADERS
            FILES
                collision_configuration.hpp
)

target_link_libraries("geometry_collisions_vparmu"
    INTERFACE
        gslx::collisions
)

add_library("gslx::geometry_collisions_vparmu" ALIAS "geometry_collisions_vparmu")

install(TARGETS geometry_collisions_vparmu EXPORT GyselalibxxTargets FILE_SET HEADERS)
