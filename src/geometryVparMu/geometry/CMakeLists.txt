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
    gyselalibxx::speciesinfo
    gyselalibxx::utils
)
add_library("gyselalibxx::geometry_vparmu" ALIAS "geometry_vparmu")

install(TARGETS geometry_vparmu EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryVparMu/geometry)
