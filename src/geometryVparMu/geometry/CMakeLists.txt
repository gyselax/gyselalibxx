# SPDX-License-Identifier: MIT

add_library("geometry_vparmu" INTERFACE)
target_sources("geometry_vparmu"
    INTERFACE
        FILE_SET HEADERS
            FILES
                geometry_vpar_mu.hpp
)
target_link_libraries("geometry_vparmu" INTERFACE
    DDC::core
    gslx::speciesinfo
    gslx::utils
)
add_library("gslx::geometry_vparmu" ALIAS "geometry_vparmu")

install(TARGETS geometry_vparmu EXPORT GyselalibxxTargets FILE_SET HEADERS)
