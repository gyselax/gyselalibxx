
add_library(quadrature INTERFACE)
target_include_directories(quadrature
    INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

target_link_libraries(quadrature INTERFACE
    DDC::core
    gslx::utils
)

add_library(gslx::quadrature ALIAS quadrature)
