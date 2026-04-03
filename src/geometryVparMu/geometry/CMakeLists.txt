# SPDX-License-Identifier: MIT

add_library("geometry_vparmu" INTERFACE)

target_compile_features("geometry_vparmu" INTERFACE cxx_std_20)
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
