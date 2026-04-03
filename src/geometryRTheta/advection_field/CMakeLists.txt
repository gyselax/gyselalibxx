

add_library("advection_field_RTheta" INTERFACE)

target_compile_features("advection_field_RTheta" INTERFACE cxx_std_20)

target_sources("advection_field_RTheta"
    INTERFACE
        FILE_SET HEADERS
            FILES
                advection_field_rtheta.hpp
)

target_link_libraries("advection_field_RTheta" INTERFACE
    DDC::core
    gslx::coord_transformations
    gslx::data_types
    gslx::pde_solvers
    gslx::poisson_RTheta
    gslx::speciesinfo
    gslx::utils
)
add_library("gslx::advection_field_RTheta" ALIAS "advection_field_RTheta")

install(TARGETS advection_field_RTheta EXPORT GyselalibxxTargets FILE_SET HEADERS)
