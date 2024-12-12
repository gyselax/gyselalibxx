


add_library(timestepper INTERFACE)
target_include_directories(timestepper
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(timestepper INTERFACE
    DDC::core
    gslx::data_types
    gslx::multipatch_data_types
    gslx::multipatch_utils
    gslx::utils
)
add_library(gslx::timestepper ALIAS timestepper)
