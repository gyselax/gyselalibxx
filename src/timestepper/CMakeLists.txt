


add_library(timestepper INTERFACE)
target_include_directories(timestepper
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(timestepper INTERFACE
    DDC::DDC
    sll::splines
    vcx::utils
)
add_library(vcx::timestepper ALIAS timestepper)
