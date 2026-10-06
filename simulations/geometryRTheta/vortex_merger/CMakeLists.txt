

add_executable(vortex_merger
     vortex_merger.cpp
)
target_link_libraries(vortex_merger
    PUBLIC
        DDC::core
        DDC::pdi

        paraconf::paraconf

        PDI::PDI_C

        gyselalibxx::advection
        gyselalibxx::geometry_RTheta
        gyselalibxx::io
        gyselalibxx::interpolation
        gyselalibxx::paraconfpp
        gyselalibxx::poisson_RTheta
        gyselalibxx::quadrature
        gyselalibxx::simulation_utils
        gyselalibxx::time_integration_RTheta
        gyselalibxx::utils
        gyselalibxx::vortex_merger_init_eq
)

