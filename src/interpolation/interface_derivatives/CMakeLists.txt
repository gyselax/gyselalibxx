# SPDX-License-Identifier: MIT

add_library("interpolation_interface_derivatives" INTERFACE)

target_compile_features("interpolation_interface_derivatives" INTERFACE cxx_std_20)

target_sources("interpolation_interface_derivatives"
    INTERFACE
        FILE_SET HEADERS
            FILES
                deriv_reconstruction_coeffs.hpp
)

target_link_libraries("interpolation_interface_derivatives"
    INTERFACE
        DDC::core
        gyselalibxx::utils
)

add_library("gyselalibxx::interpolation_interface_derivatives" ALIAS "interpolation_interface_derivatives")

install(TARGETS interpolation_interface_derivatives EXPORT gyselalibxxTargets FILE_SET HEADERS)
