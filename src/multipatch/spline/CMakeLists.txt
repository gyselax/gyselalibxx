

add_library("multipatch_spline" INTERFACE)
target_include_directories("multipatch_spline"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("multipatch_spline" INTERFACE
    DDC::core
    gslx::multipatch_data_types
    gslx::multipatch_connectivity
    gslx::utils

)
add_library("gslx::multipatch_spline" ALIAS "multipatch_spline")
