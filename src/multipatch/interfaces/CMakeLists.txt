cmake_minimum_required(VERSION 3.15)


add_library("multipatch_interface" INTERFACE)
target_include_directories("multipatch_interface"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("multipatch_interface" INTERFACE
    DDC::DDC
    gslx::utils
)
add_library("gslx::multipatch_interface" ALIAS "multipatch_interface")
