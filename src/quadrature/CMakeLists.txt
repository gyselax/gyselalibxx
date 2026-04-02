
add_library(quadrature STATIC
    gauss_legendre_integration.cpp
    )

target_sources("quadrature"
    INTERFACE
        FILE_SET HEADERS
            FILES
                gauss_legendre_integration.hpp
                neumann_spline_quadrature.hpp
                quadrature_coeffs_nd.hpp
                quadrature.hpp
                simpson_quadrature.hpp
                spline_quadrature.hpp
                trapezoid_quadrature.hpp
                volume_quadrature_nd.hpp
)

target_link_libraries(quadrature PUBLIC
    DDC::core
    gslx::coord_transformations
    gslx::utils
)

add_library(gslx::quadrature ALIAS quadrature)

install(TARGETS quadrature EXPORT GyselalibxxTargets FILE_SET HEADERS)
