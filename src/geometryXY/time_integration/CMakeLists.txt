# SPDX-License-Identifier: MIT

add_library("predcorr_rk2_xy" INTERFACE)

target_sources("predcorr_rk2_xy"
    INTERFACE
        FILE_SET HEADERS
            FILES
                predcorr_RK2.hpp
)

target_link_libraries("predcorr_rk2_xy" INTERFACE
    DDC::core
    gslx::advection
    gslx::interpolation
    gslx::geometry_XY
    gslx::pde_solvers
    gslx::timestepper
    gslx::utils
)
add_library("gslx::predcorr_rk2_xy" ALIAS "predcorr_rk2_xy")

install(TARGETS predcorr_rk2_xy EXPORT GyselalibxxTargets FILE_SET HEADERS)
