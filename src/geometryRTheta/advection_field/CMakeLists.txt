

add_library("advection_field_RTheta" INTERFACE)
target_include_directories("advection_field_RTheta"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("advection_field_RTheta" INTERFACE
    DDC::DDC
    sll::SLL
    gslx::data_types
    gslx::speciesinfo
    gslx::utils
    gslx::poisson_RTheta
)
add_library("gslx::advection_field_RTheta" ALIAS "advection_field_RTheta")
