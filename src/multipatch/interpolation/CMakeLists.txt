# SPDX-License-Identifier: MIT

add_library("multipatch_interpolation" INTERFACE)

target_compile_features("multipatch_interpolation" INTERFACE cxx_std_20)

target_sources("multipatch_interpolation"
    INTERFACE
        FILE_SET HEADERS
            FILES
                interpolation_builder_2d_collection.hpp
                interpolation_builder_collection.hpp
)

target_link_libraries("multipatch_interpolation"
    INTERFACE
        DDC::core
        gyselalibxx::interpolation
        gyselalibxx::multipatch_connectivity
        gyselalibxx::multipatch_data_types
        gyselalibxx::utils
)

add_library("gyselalibxx::multipatch_interpolation" ALIAS "multipatch_interpolation")

install(TARGETS multipatch_interpolation EXPORT gyselalibxxTargets FILE_SET HEADERS)
