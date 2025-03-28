
add_library("multipatch_data_types" INTERFACE)
target_include_directories("multipatch_data_types"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("multipatch_data_types" INTERFACE
    DDC::core
    DDC::splines

    gslx::data_types
    gslx::utils
)
add_library("gslx::multipatch_data_types" ALIAS "multipatch_data_types")
