# SPDX-License-Identifier: MIT

add_library("geometry_XY" INTERFACE)

target_compile_features("geometry_XY" INTERFACE cxx_std_20)

target_sources("geometry_XY"
    INTERFACE
        FILE_SET HEADERS
            FILES
                geometry_xy.hpp
)

target_link_libraries("geometry_XY" INTERFACE
    DDC::core
    DDC::splines
    gyselalibxx::utils
)
add_library("gyselalibxx::geometry_XY" ALIAS "geometry_XY")

install(TARGETS geometry_XY EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXY/geometry)
