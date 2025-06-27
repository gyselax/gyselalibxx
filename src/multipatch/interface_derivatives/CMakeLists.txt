
add_library("multipatch_interface_derivatives" INTERFACE)
target_include_directories("multipatch_interface_derivatives"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("multipatch_interface_derivatives" INTERFACE
    DDC::core
    gslx::multipatch_data_types
    gslx::multipatch_connectivity
    gslx::utils
)
add_library("gslx::multipatch_interface_derivatives" ALIAS "multipatch_interface_derivatives")
