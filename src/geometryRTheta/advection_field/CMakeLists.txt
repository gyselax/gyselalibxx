

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
    gyselalibxx::coord_transformations
    gyselalibxx::data_types
    gyselalibxx::pde_solvers
    gyselalibxx::poisson_RTheta
    gyselalibxx::speciesinfo
    gyselalibxx::utils
)
add_library("gyselalibxx::advection_field_RTheta" ALIAS "advection_field_RTheta")

install(TARGETS advection_field_RTheta EXPORT GyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryRTheta/advection_field)
