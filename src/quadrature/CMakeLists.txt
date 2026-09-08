
add_library(quadrature STATIC
    gauss_legendre_integration.cpp
    )

target_compile_features("quadrature" PUBLIC cxx_std_20)

target_sources("quadrature"
    INTERFACE
        FILE_SET HEADERS
            FILES
                gauss_legendre_integration.hpp
                neumann_spline_quadrature.hpp
                quadrature.hpp
                quadrature_coeffs_nd.hpp
                simpson_quadrature.hpp
                spline_quadrature.hpp
                trapezoid_quadrature.hpp
                volume_quadrature_nd.hpp
)

target_link_libraries(quadrature PUBLIC
    DDC::core
    gyselalibxx::coord_transformations
    gyselalibxx::utils
)

add_library(gyselalibxx::quadrature ALIAS quadrature)

install(TARGETS quadrature EXPORT gyselalibxxTargets FILE_SET HEADERS)
