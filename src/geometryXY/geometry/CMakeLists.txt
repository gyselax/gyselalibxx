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
    gslx::utils
)
add_library("gslx::geometry_XY" ALIAS "geometry_XY")

install(TARGETS geometry_XY EXPORT GyselalibxxTargets FILE_SET HEADERS)
