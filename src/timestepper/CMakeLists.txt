


add_library(timestepper INTERFACE)
target_include_directories(timestepper
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(timestepper INTERFACE
    DDC::DDC
    sll::splines
    gslx::utils
)
add_library(gslx::timestepper ALIAS timestepper)
