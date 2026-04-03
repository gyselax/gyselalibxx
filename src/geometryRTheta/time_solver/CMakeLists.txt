# SPDX-License-Identifier: MIT


add_library(time_integration_RTheta INTERFACE)

target_compile_features("time_integration_RTheta" INTERFACE cxx_std_20)

target_sources("time_integration_RTheta"
    INTERFACE
        FILE_SET HEADERS
            FILES
                bsl_predcorr.hpp
                bsl_predcorr_second_order_explicit.hpp
                bsl_predcorr_second_order_implicit.hpp
                itimesolver.hpp
)

target_link_libraries(time_integration_RTheta
    INTERFACE
        DDC::core
        gslx::advection
        gslx::advection_field_RTheta
        gslx::geometry_RTheta
        gslx::poisson_RTheta
        gslx::utils
        gslx::timestepper
)
add_library(gslx::time_integration_RTheta ALIAS time_integration_RTheta)

install(TARGETS time_integration_RTheta EXPORT GyselalibxxTargets FILE_SET HEADERS)
