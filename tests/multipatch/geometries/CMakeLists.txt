
add_library("multipatch_geometries" INTERFACE)
target_include_directories("multipatch_geometries"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("multipatch_geometries" INTERFACE
    DDC::core
    gyselalibxx::multipatch_connectivity
    gyselalibxx::utils

)
add_library("gyselalibxx::multipatch_geometries" ALIAS "multipatch_geometries")
