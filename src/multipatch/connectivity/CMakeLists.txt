

add_library("multipatch_connectivity" INTERFACE)
target_include_directories("multipatch_connectivity"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("multipatch_connectivity" INTERFACE
    DDC::DDC
    gslx::utils
    gslx::data_types
)
add_library("gslx::multipatch_connectivity" ALIAS "multipatch_connectivity")
