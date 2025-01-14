
add_library(quadrature STATIC
    gauss_legendre_integration.cpp
    )
target_include_directories(quadrature
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

target_link_libraries(quadrature PUBLIC
    DDC::core
    sll::SLL
    gslx::utils
)

add_library(gslx::quadrature ALIAS quadrature)
