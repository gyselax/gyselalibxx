

add_executable(vortex_merger
     vortex_merger.cpp
)
target_link_libraries(vortex_merger
    PUBLIC
        DDC::core
        DDC::pdi

        paraconf::paraconf

        PDI::pdi

        gslx::advection
        gslx::geometry_RTheta
        gslx::io
        gslx::interpolation
        gslx::paraconfpp
        gslx::poisson_RTheta
        gslx::quadrature
        gslx::simulation_utils
        gslx::time_integration_RTheta
        gslx::utils
        gslx::vortex_merger_init_eq
)

