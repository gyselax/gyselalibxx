# SPDX-License-Identifier: MIT

add_library("predcorr_rk2_xy" INTERFACE)

target_compile_features("predcorr_rk2_xy" INTERFACE cxx_std_20)

target_sources("predcorr_rk2_xy"
    INTERFACE
        FILE_SET HEADERS
            FILES
                predcorr_RK2.hpp
)

target_link_libraries("predcorr_rk2_xy" INTERFACE
    DDC::core
    gyselalibxx::advection
    gyselalibxx::interpolation
    gyselalibxx::geometry_XY
    gyselalibxx::pde_solvers
    gyselalibxx::timestepper
    gyselalibxx::utils
)
add_library("gyselalibxx::predcorr_rk2_xy" ALIAS "predcorr_rk2_xy")

install(TARGETS predcorr_rk2_xy EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXY/time_integration)
