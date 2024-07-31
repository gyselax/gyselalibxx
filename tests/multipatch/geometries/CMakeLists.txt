
add_library("multipatch_geometries" INTERFACE)
target_include_directories("multipatch_geometries"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("multipatch_geometries" INTERFACE
    DDC::DDC
    gslx::multipatch_connectivity
    gslx::utils

)
add_library("gslx::multipatch_geometries" ALIAS "multipatch_geometries")
